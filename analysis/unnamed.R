library(dplyr)
library(stringr)


stats <- read.table("results/LMD.sigDE_stats.txt", sep = "\t",
                    header = TRUE, stringsAsFactors = FALSE)

str(stats)
# 'data.frame':	2283 obs. of  4 variables:

upReg <- stats %>% filter(logFC > 0)

upRegBoth <- upReg[duplicated(upReg$cds), "cds"]

upReg[upReg$cds %in% upRegBoth, "focus_term"] <- "Both"

colnames(upReg)[4] <- "Up_Regulation"


blast <- read.table("data/SCB.upRegDSC.17Jun.blastpTAIR10.txt",
                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)

blast <- blast %>% select(qseqid, sseqid, evalue, salltitles)

blast$sseqid <- blast$sseqid %>% 
  str_replace(".\\d$", "")

colnames(blast)[2] <- "TAIR_ID"

str(blast)
# 'data.frame':	1108 obs. of  4 variables:


table1 <- right_join(upReg, blast, by = c("cds" = "qseqid"))

str(table1)
# 'data.frame':	1171 obs. of  7 variables:

write.table(table1, "results/table1.txt", 
            quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)


secondary_cellwall_tairIDs <- 
  scan("results/secondary_cellwall_tairIDs.txt", what = "character")

table2 <- table1 %>% filter(TAIR_ID %in% secondary_cellwall_tairIDs)

str(table2)
# 'data.frame':	79 obs. of  7 variables:

write.table(table2, "results/table2.txt", 
            quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)
