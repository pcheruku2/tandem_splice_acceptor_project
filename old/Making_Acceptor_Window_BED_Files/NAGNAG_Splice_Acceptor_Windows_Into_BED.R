rm(list = ls())
library(dplyr)
library(magrittr)
library(tidyr)
library(data.table)
library(pbmcapply)

homeDir <- "/Users/frostfg/Desktop/Files-and-such/NAGNAG/Genome-wide_NAGNAG"
setwd(homeDir)

AcceptorWindows <- fread("KnownGenes_UCSC_GRCh37_MANE_Transcripts_Splice_Acceptor_Windows.tsv")
colnames(AcceptorWindows)[1] <- "refseqID"
  
bedFormat <- 
  function(id,sourceData) {
    # id <- "NM_032785.4"
    # sourceData <- AcceptorWindows
    starts <- unlist(strsplit(sourceData$ExonWindowStart[sourceData$refseqID == id],","))
    ends <- unlist(strsplit(sourceData$ExonWindowEnd[sourceData$refseqID == id],","))
    outExon <- 
      data.frame("chrom" = sourceData$chrom[sourceData$refseqID == id],
                 "chromStart" = as.double(starts) - 1, # Convert to BED coordinate system
                 "chromEnd" = ends,
                 "name" = sourceData$refseqID[sourceData$refseqID == id],
                 "score" = 0,
                 "strand" = sourceData$strand[sourceData$refseqID == id]
                 )
    
    starts <- unlist(strsplit(sourceData$IntronWindowStart[sourceData$refseqID == id],","))
    ends <- unlist(strsplit(sourceData$IntronWindowEnd[sourceData$refseqID == id],","))
    outIntron <- 
      data.frame("chrom" = sourceData$chrom[sourceData$refseqID == id],
                 "chromStart" = as.double(starts) - 1, # Convert to BED coordinate system
                 "chromEnd" = ends,
                 "name" = sourceData$refseqID[sourceData$refseqID == id],
                 "score" = 0,
                 "strand" = sourceData$strand[sourceData$refseqID == id]
                 )
    out <- 
      list("ExonWindow" = outExon,
           "IntronWindow" = outIntron
           )
    rm(outExon)
    rm(outIntron)
    return(out)
    }
 
bedWindows <-
  AcceptorWindows$refseqID %>%
  unique %>%
  as.list %>%
  pbmclapply(bedFormat,
             sourceData = AcceptorWindows,
             mc.cores = detectCores() - 1
             )

bedIntronWindows <-
  bedWindows %>% 
  lapply(`[[`,"IntronWindow") %>% 
  bind_rows
bedExonWindows <-
  bedWindows %>% 
  lapply(`[[`,"ExonWindow") %>% 
  bind_rows   

fwrite(bedIntronWindows,
       "GRCh37_MANE_Transcripts_Splice_Acceptor_Intronic_Windows.bed",
       sep = "\t",
       col.names = FALSE
       )   
fwrite(bedExonWindows,
       "GRCh37_MANE_Transcripts_Splice_Acceptor_Exonic_Windows.bed",
       sep = "\t",
       col.names = FALSE
       )   
