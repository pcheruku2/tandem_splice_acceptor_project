#!/bin/R
library(tidyverse)
library(Biostrings)
library(pbmcapply)
library(magrittr)

homeDir <- "/path/to/output"
chrDataDir <- "/path/to/chromosome/fastA/files"

setwd(chrDataDir)

intronFastaFiles <- list.files()[grep("Intronic",list.files())]

allPolyPyrTracts <-
	intronFastaFiles %>%
	as.list %>%
	pbmclapply(function(file){
			# file <- "chrY_-_Intronic_regions.fa"
			fastA <- readDNAStringSet(file)
			fastA <- as.character(fastA)
			fastApyr <- str_replace_all(fastA,"[CT]","1")
			fastApyr <- str_replace_all(fastApyr,"[AG]","0")
			
			pPyrLen <- str_extract_all(fastApyr, "1+")
			pyrPatterns <- 
				lapply(pPyrLen,
				       function(x){
					       out <- x[which.max(nchar(x))]
					       out <- ifelse(length(out) == 0,NA,out)
					       return(out)
				       }
				       )
			pPyrLen <- 
				lapply(pPyrLen,
				       function(x){max(nchar(x))}
				       )
			pPyrLen <- unlist(pPyrLen)
			pPyrLen[is.infinite(pPyrLen)] <- NA
			
			pPyrLoc <- 
				str_locate_all(fastApyr,
					       unlist(pyrPatterns)
					       ) %>%
				lapply(function(x){
					       if(nrow(x) > 1){
						       x <- x[x[,"end"] == max(x[,"end"]),]
					       } else {}
					       x <- x - 33
					       return(x)
					       }
				       )
			pPyrLoc <- do.call(rbind,pPyrLoc)
			out <- 
				data.frame("Sequence" = fastA,
					   "seqName" = names(fastA),
					   "Pyrimidines" = fastApyr,
					   "PolyPyrTract_Length" = pPyrLen,
					   pPyrLoc
					   )
			},
		   mc.cores = 24
		   )

allPolyPyrTracts %<>% bind_rows
colnames(allPolyPyrTracts)[colnames(allPolyPyrTracts) %in% c("start","end")] <- c("TractStart","TractEnd")
rownames(allPolyPyrTracts) <- NULL
allPolyPyrTracts %<>%
	mutate(Chromosome = str_remove(str_extract(seqName,"^.*:"),":"),
	       WindowStart = str_remove(str_extract(seqName,"[0-9]*-"),"-"),
	       WindowEnd = str_remove(str_extract(seqName,"-[0-9]*"),"-")
	       ) %>%
	select(Chromosome,WindowStart,WindowEnd,
	       PolyPyrTract_Length,TractStart,TractEnd,
	       Sequence,Pyrimidines
	       )

allPolyPyrTracts %>%
	write_tsv("GRCh37_MANE_Transcripts_PolyPyrimidine_Tract_Data.tsv")

