rm(list = ls())
library(dplyr)
library(magrittr)
library(stringr)
library(data.table)
library(pbmcapply)

homeDir <- "/path/to/splice/acceptor/region/files"
setwd(homeDir)

ManeTranscripts <- 
  fread("KnownGenes_UCSC_GRCh37_MANE_Transcripts_Only.tsv") %>%
  filter(exonCount > 1) %>%
  select(name,GeneSymbol,chrom,strand,
         txStart,txEnd,cdsStart,cdsEnd,
         exonCount,exonStarts,exonEnds
         ) 

WindowFinder <- 
  function(transcript, genedata) {
    # Test input, source data: genedata <- ManeTranscripts
    # Test input, + strand: transcript <- "NM_052998.4"
    # Test input, - strand: transcript <- "NM_014288.5"
    geneInfo <- genedata[genedata$name == transcript,]

    if(geneInfo$strand == "+") {
      # Normal + strand direction
      # Start by defining the coding exons
      exonPositionsStart <- as.double(unlist(str_split(gsub('.{1}$','',geneInfo$exonStarts),",")))
      # This is the first exon that fully contains the cds, i.e. the first exon with a splice acceptor site
      firstCDSexon <- 
        match(head(exonPositionsStart[exonPositionsStart > geneInfo$cdsStart],
                   n = 1
                   ),
              exonPositionsStart
              )
      
      # This is the last exon that fully contains the cds, i.e. the last exon with a splice acceptor site
      lastCDSexon <- 
        match(tail(exonPositionsStart[exonPositionsStart < geneInfo$cdsEnd],
                   n = 1
                   ),
              exonPositionsStart
              )
      if(length(firstCDSexon) == 0 | length(lastCDSexon) == 0) {
        out <- NULL
        windowStart <- "none"
        windowEnd <- "none"
      } else {
        IntronWindowStart <- exonPositionsStart[firstCDSexon:lastCDSexon] - 31
        IntronWindowEnd <- exonPositionsStart[firstCDSexon:lastCDSexon] + 1
        
        ExonWindowStart <- exonPositionsStart[firstCDSexon:lastCDSexon] + 0
        ExonWindowEnd <- exonPositionsStart[firstCDSexon:lastCDSexon] + 30
        
        
        out <-
          data.frame(geneInfo,
                     "IntronWindowStart" = paste(as.character(IntronWindowStart), collapse = ","),
                     "IntronWindowEnd" = paste(as.character(IntronWindowEnd), collapse = ","),
                     "ExonWindowStart" = paste(as.character(ExonWindowStart), collapse = ","),
                     "ExonWindowEnd" = paste(as.character(ExonWindowEnd), collapse = ",")
                     )
      }
      
      rm(list = c("exonPositionsStart","firstCDSexon","lastCDSexon","IntronWindowStart","IntronWindowEnd","ExonWindowStart","ExonWindowEnd"))
    } else {
      # I actually have to think for this one
      # Start by defining the coding exons
      exonPositionsEnd <- as.double(unlist(str_split(gsub('.{1}$','',geneInfo$exonEnds),",")))
      
      # This is the first exon that fully contains the cds, i.e. the first exon with a splice acceptor site
      firstCDSexon <- 
        match(tail(exonPositionsEnd[exonPositionsEnd < geneInfo$cdsEnd],
                   n = 1
                   ),
              exonPositionsEnd
        )
      # This is the last exon that fully contains the cds, i.e. the last exon with a splice acceptor site
      lastCDSexon <- 
        match(head(exonPositionsEnd[geneInfo$cdsStart < exonPositionsEnd],
                   n = 1
                   ),
              exonPositionsEnd
              )
      if(length(firstCDSexon) == 0 | length(lastCDSexon) == 0) {
        out <- NULL
        windowStart <- "none"
        windowEnd <- "none"
        
      } else {
        IntronWindowStart <- exonPositionsEnd[firstCDSexon:lastCDSexon] + 32
        IntronWindowEnd <- exonPositionsEnd[firstCDSexon:lastCDSexon] + 0
        
        ExonWindowStart <- exonPositionsEnd[firstCDSexon:lastCDSexon] + 1
        ExonWindowEnd <- exonPositionsEnd[firstCDSexon:lastCDSexon] - 29
        out <-
          # This Shit is going to be backwards for the - strand so that the bed file isn't jacked up
          data.frame(geneInfo,
                     "IntronWindowStart" =  paste(as.character(IntronWindowEnd), collapse = ","),
                     "IntronWindowEnd" = paste(as.character(IntronWindowStart), collapse = ","),
                     "ExonWindowStart" = paste(as.character(ExonWindowEnd), collapse = ","),
                     "ExonWindowEnd" = paste(as.character(ExonWindowStart), collapse = ",")
                     )
      }
      
      rm(list = c("exonPositionsStart","firstCDSexon","lastCDSexon","IntronWindowStart","IntronWindowEnd","ExonWindowStart","ExonWindowEnd"))
    }
    return(out)
  }


WindowData <-
  ManeTranscripts$name %>%
  as.list %>%
  pbmclapply(WindowFinder,
         genedata = ManeTranscripts,
         mc.cores = detectCores() -1
         )

WindowData %>% 
  bind_rows %>%
  select(name,GeneSymbol,chrom,strand,
         IntronWindowStart,IntronWindowEnd,
         ExonWindowStart,ExonWindowEnd
         ) %>% 
  fwrite("KnownGenes_UCSC_GRCh37_MANE_Transcripts_Splice_Acceptor_Windows.tsv",
         sep = "\t"
         )
