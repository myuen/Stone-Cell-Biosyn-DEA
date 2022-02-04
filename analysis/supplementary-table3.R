library(dplyr)
library(stringr)


# R script to write out annotation for supplementary tables

# Read DEA results
stats <- read.table("results/LMD.sigDE_stats.txt", sep = "\t",
                    header = TRUE, stringsAsFactors = FALSE)

str(stats)
# 'data.frame':	2283 obs. of  4 variables:


# Read tab-delimited BLAST output
blast <- read.delim("data/SCB.upRegDSC.17Jun.blastpNR.txt",
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)

blast <- blast %>% 
  select(qseqid, sseqid, evalue, pident, ppos, salltitles)


# Join up-regulated DEA statistics with BLAST annotations
tableS3 <- left_join(stats, blast, by = c("cds" = "qseqid"))

str(tableS3)
# 'data.frame':	13390 obs. of  9 variables:


tableS3 %>% group_by(focus_term, logFC > 0) %>% 
  summarise(count = n())
#   focus_term `logFC > 0` count
# 1 R_cType    FALSE         495
# 2 R_cType    TRUE        10347
# 3 S_cType    FALSE         451
# 4 S_cType    TRUE         2097


tableS3 %>%
  filter(logFC > 0) %>% 
  group_by(focus_term, .drop = FALSE) %>% 
  group_walk(.keep = TRUE, ~ write.table(.x, file = paste0("results/", .y, ".upReg.txt"),
                           quote = FALSE, row.names = FALSE, sep = "\t"))

tableS3 %>%
  filter(logFC < 0) %>% 
  group_by(focus_term, .drop = FALSE) %>% 
  group_walk(.keep = TRUE, ~ write.table(.x, file = paste0("results/", .y, ".downReg.txt"),
                                         quote = FALSE, row.names = FALSE, sep = "\t"))
