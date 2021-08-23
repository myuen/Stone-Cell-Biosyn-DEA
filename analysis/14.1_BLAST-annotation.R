library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)


### Read BLAST results for all DE contigs
blast <- read.table("results//SCB.upRegDSC.17Jun.blastpNR.txt", 
                    header = TRUE, sep = "\t", 
                    quote = "", stringsAsFactors = FALSE)

blast.tophit <- blast %>% nest(data = -qseqid)

blast.tophit$data <- blast.tophit$data %>%
  map(function(x){
  x[1,]})

blast.tophit <- blast.tophit %>% unnest(-qseqid)

str(blast.tophit)
# tibble [1,161 × 21] (S3: tbl_df/tbl/data.frame)


sctc.s_peakedT2 <- 
  scan("results/SCTC.S_PeakedT2.txt", what = "character")
sctc.s_peakedT3 <- 
  scan("results/SCTC.S_PeakedT3.txt", what = "character")
sctc.r_peakedT2 <- 
  scan("results/SCTC.R_PeakedT2.txt", what = "character")
sctc.r_peakedT3 <- 
  scan("results/SCTC.R_PeakedT3.txt", what = "character")



### S T2
sctc.s_peakedT2.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.s_peakedT2)

str(sctc.s_peakedT2.tophit)
# tibble [25 × 21] (S3: tbl_df/tbl/data.frame)

write.table(sctc.s_peakedT2.tophit, 
            "results/SCTC.S_PeakedT2.blastTopHit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### S T3
sctc.s_peakedT3.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.s_peakedT3)

str(sctc.s_peakedT3.tophit)
# tibble [24 × 21] (S3: tbl_df/tbl/data.frame)

write.table(sctc.s_peakedT3.tophit, 
            "results/SCTC.S_PeakedT3.blastTopHit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T2
sctc.r_peakedT2.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.r_peakedT2)

str(sctc.r_peakedT2.tophit)
# tibble [143 × 21] (S3: tbl_df/tbl/data.frame)

write.table(sctc.r_peakedT2.tophit, 
            "results/SCTC.R_PeakedT2.blastTopHit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T3
sctc.r_peakedT3.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.r_peakedT3)

write.table(sctc.r_peakedT3.tophit, 
            "results/SCTC.R_PeakedT3.blastTopHit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
