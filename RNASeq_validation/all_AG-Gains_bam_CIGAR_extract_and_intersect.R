#!/bin/R

args <- commandArgs(trailingOnly=TRUE)
variantBamFile <- args[1]
individualCaseID <- args[2]

suppressMessages(library(stringr))
suppressMessages(library(readr,include.only = c("read_tsv","cols","write_tsv")))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr,include.only = c("%>%","%<>%")))
suppressMessages(library(stringr))
suppressMessages(library(GenomicAlignments,include.only = c("readGAlignments","njunc","cigar")))
suppressMessages(library(GenomicRanges))
suppressWarnings(suppressMessages(library(sqldf)))


# Found some code to map the junctions, seems to work well!
# Note: the comments are lightly modified from the original code for clarity
# http://malcook-gedanken.blogspot.com/2011/06/tabulating-reads_09.html

bamObject <- readGAlignments(variantBamFile)
bamObjectJunctions <- bamObject[njunc(bamObject) > 0]
bamObjectSoftClip <- bamObject[grepl("S",cigar(bamObject))]
bamObjectSoftClip <- as.data.frame(bamObjectSoftClip)
rm(bamObject)


# Makes a GRanges, with one component per alignment (which looses spliced internals).
bamObjectJunctions.grg <- granges(bamObjectJunctions)

# Makes a GRangesList of GRanges, one per alignment
bamObjectJunctions.grglist <- grglist(bamObjectJunctions)

# Extracts the gapped regions (presumably intronic), as GRangesList, indexed by read.
J2J.grglist <- psetdiff(bamObjectJunctions.grg,bamObjectJunctions.grglist)

# Extract the gapped regions, as genome wide GRanges
J2J.GRanges <- unlist(J2J.grglist)

# Offset positions to include the last base of upstream exon and first base of downstream exon
start(J2J.GRanges) <- start(J2J.GRanges) - 1
end(J2J.GRanges)   <- end(J2J.GRanges)   + 1

# Coerce to data.frame so we can query using sqldf.
J2JDF <- as.data.frame(J2J.GRanges)
# Count the number of reads grouped by chromosome,start,end,strand
tabulatedGaps <-
	paste(
		"SELECT seqnames AS chr,start,end, count(*) AS n",
                "FROM J2JDF",
                "GROUP BY seqnames,start,end"
		) %>%
	sqldf

# Collect variant information from input file name
caseBamFile <-
	variantBamFile %>% 
	str_remove("/path/to/AG-gain/BAM/files/") %>%
        str_remove("_AG-gain-supporting_reads\\/chr.*$")

variantBamNoPath <-
        variantBamFile %>%  
        str_remove(".*supporting_reads\\/") 
variantInfo <-
	variantBamNoPath %>%
	str_remove("_AG-Gain.*$")
variantChr <- 
	variantInfo %>%
	str_extract("chr.{1,2}_") %>%
	str_remove("_")
variantPos <-
        variantInfo %>%
        str_extract("_[0-9]*_") %>%
        str_remove_all("_")

bamRangeEnd <-
       variantBamNoPath %>%
       str_extract("[0-9]{1,20}-[0-9]{1,20}") %>%
       str_remove("[0-9]{1,20}-") %>%
       as.double
bamRangeStart <- 
	       variantBamNoPath %>%     
	       str_extract("[0-9]{1,20}-[0-9]{1,20}") %>%
	       str_remove("-[0-9]{1,20}") %>%
	       as.double

outDir <- 
	paste0(caseBamFile,"_AG-gain-supporting_reads_CIGAR_Out") %>%
        str_replace("AG-Gain_BAM_Files","all_AG-Gain_CIGAR_Decoding")


# Load Reference NAGNAG DF for this case and variant and mark the genomic position of the new exon start
refNagnagDF <-
	"/path/to/All_NAGNAG_Vars_With_NewExonStart.tsv" %>%
	read_tsv(col_types = cols(),
		 num_threads = 1
		 ) %>%
	filter(
		individualID == individualCaseID,
		Chromosome == variantChr,
		Position == variantPos
		) %>%
	mutate(
                refExonGenomicStart = NA,
		newExonGenomicStart = NA,
		junctionNewExonReads = 0,
		junctionRefExonReads = 0,
		softClipNewExonReads = NA
                                        #Strand == "+",
                                        #nrow(bamObjectSoftClip[bamObjectSoftClip$start == newExonGenomicStart,]),
                                        #nrow(bamObjectSoftClip[bamObjectSoftClip$end == newExonGenomicStart,])
		) %>%
	select(
		Chromosome,Position,Strand,Ref,Alt,
		ORF_Ref,ORF_Alt,RefTrimer,ORF_AltTrimer,
		VariantType,AGPositionWithinTrimer,BoundaryPosition,newExonPosition,
		refExonGenomicStart,newExonGenomicStart,
		junctionRefExonReads,junctionNewExonReads,softClipNewExonReads,
		individualID
		) %>%
	as.data.frame
if(nrow(refNagnagDF) > 1) {
        print("More than one variant listed, taking only the first row:")
        print(refNagnagDF)
        refNagnagDF <- refNagnagDF[1,]
}
# Set what the ref/alt exon boundaries should be based on intronic/exonic and +/- strand
# + strand exonic variants - Checked
if(refNagnagDF$Strand == "+" & refNagnagDF$BoundaryPosition > 0) {
	refNagnagDF$refExonGenomicStart <- bamRangeStart
	refNagnagDF$newExonGenomicStart <- bamRangeEnd
	refNagnagDF$softClipNewExonReads <-
		ifelse(refNagnagDF$Strand == "+",
                       nrow(bamObjectSoftClip[bamObjectSoftClip$start == refNagnagDF$newExonGenomicStart,]),
                       nrow(bamObjectSoftClip[bamObjectSoftClip$end == refNagnagDF$newExonGenomicStart,])
		       )
	
	refNagnagDF$junctionRefExonReads <- sum(tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$refExonGenomicStart])
        if(length(tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$newExonGenomicStart]) > 0){
                refNagnagDF$junctionNewExonReads <- (tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$newExonGenomicStart])
        }
}

# + strand intronic variants - Checked
if(refNagnagDF$Strand == "+" & refNagnagDF$BoundaryPosition < 0) {
        refNagnagDF$refExonGenomicStart <- bamRangeEnd
        refNagnagDF$newExonGenomicStart <- bamRangeStart + 1
	refNagnagDF$softClipNewExonReads <-
                ifelse(refNagnagDF$Strand == "+",
                       nrow(bamObjectSoftClip[bamObjectSoftClip$start == refNagnagDF$newExonGenomicStart,]),
                       nrow(bamObjectSoftClip[bamObjectSoftClip$end == refNagnagDF$newExonGenomicStart,])
                       )

	refNagnagDF$junctionRefExonReads <- sum(tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$refExonGenomicStart])
        if(length(tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$newExonGenomicStart]) > 0){
                refNagnagDF$junctionNewExonReads <- (tabulatedGaps$n[tabulatedGaps$end == refNagnagDF$newExonGenomicStart])
        }
}

# - strand exonic variants - Checked
if(refNagnagDF$Strand == "-" & refNagnagDF$BoundaryPosition > 0) {
        refNagnagDF$refExonGenomicStart <- bamRangeEnd
	refNagnagDF$newExonGenomicStart <- bamRangeStart
	refNagnagDF$softClipNewExonReads <-
                ifelse(refNagnagDF$Strand == "+",
                       nrow(bamObjectSoftClip[bamObjectSoftClip$start == refNagnagDF$newExonGenomicStart,]),
                       nrow(bamObjectSoftClip[bamObjectSoftClip$end == refNagnagDF$newExonGenomicStart,])
                       )

	refNagnagDF$junctionRefExonReads <- sum(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$refExonGenomicStart])
	if(length(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$newExonGenomicStart]) > 0){
		refNagnagDF$junctionNewExonReads <- sum(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$newExonGenomicStart])
	}
}
# - strand intronic variants - 
if(refNagnagDF$Strand == "-" & refNagnagDF$BoundaryPosition < 0) {
        refNagnagDF$refExonGenomicStart <- bamRangeStart
        refNagnagDF$newExonGenomicStart <- bamRangeEnd
	refNagnagDF$softClipNewExonReads <-
                ifelse(refNagnagDF$Strand == "+",
                       nrow(bamObjectSoftClip[bamObjectSoftClip$start == refNagnagDF$newExonGenomicStart,]),
                       nrow(bamObjectSoftClip[bamObjectSoftClip$end == refNagnagDF$newExonGenomicStart,])
                       )
	
	refNagnagDF$junctionRefExonReads <- sum(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$refExonGenomicStart])
        if(length(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$newExonGenomicStart]) > 0){
                refNagnagDF$junctionNewExonReads <- sum(tabulatedGaps$n[tabulatedGaps$start == refNagnagDF$newExonGenomicStart])
        }

}


system(paste("mkdir -p",outDir))
setwd(outDir)

outFileName <- paste0(individualCaseID,"_",variantInfo,"_Intronic_AG_Gain_Variant_Usage.tsv")

refNagnagDF %>% write_tsv(outFileName)

paste("There are", refNagnagDF$junctionRefExonReads,"reference and",
      refNagnagDF$junctionNewExonReads,"alt exon reads, and",
      refNagnagDF$softClipNewExonReads, "soft-clipped alt exon reads"
      ) %>%
      print




