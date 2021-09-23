library(dplyr)
library(tibble)


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb_cpm <- read.delim("results/SCB.tmm_normalized_cpm.txt")

scb_cpm <- scb_cpm %>%
  rownames_to_column("cds")

str(scb_cpm)
# 'data.frame':	26787 obs. of  13 variables:

# This is cumbersome but cannot get mutate to work 
# with rowMeans(select(starts_with(col)))
scb_cpm <- scb_cpm %>% 
  rowwise() %>% 
  mutate(R_CP = mean(c(R_CP_rep1, R_CP_rep2, R_CP_rep3)),
         R_DSC = mean(c(R_DSC_rep1, R_DSC_rep2, R_DSC_rep3)),
         S_CP = mean(c(S_CP_rep1, S_CP_rep2, S_CP_rep3)),
         S_DSC = mean(c(S_DSC_rep1, S_DSC_rep2, S_DSC_rep3)),
         .keep = "unused")

str(scb_cpm)
# rowwise_df [26,787 × 5] (S3: rowwise_df/tbl_df/tbl/data.frame)

write.table(
  scb_cpm,
  "results/SCB.rowmean_cpm.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)



# Read TMM normalized CPM from time-course (TC) experiment mapped to 
# differentially expressed CDS in LMD experiment
sctc_cpm <- read.delim("results/SCTC.tmm_normalized_cpm.txt")

sctc_cpm <- sctc_cpm %>%
  rownames_to_column("cds")

str(sctc_cpm)
# 'data.frame':	1293 obs. of  25 variables:

sctc_cpm <- sctc_cpm %>%
  rowwise() %>%
  mutate(R_T1 = mean(c(R_T1_rep1, R_T1_rep2, R_T1_rep3)),
         R_T2 = mean(c(R_T2_rep1, R_T2_rep2, R_T2_rep3)),
         R_T3 = mean(c(R_T3_rep1, R_T3_rep2, R_T3_rep3)),
         R_T4 = mean(c(R_T4_rep1, R_T4_rep2, R_T4_rep3)),
         S_T1 = mean(c(S_T1_rep1, S_T1_rep2, S_T1_rep3)),
         S_T2 = mean(c(S_T2_rep1, S_T2_rep2, S_T2_rep3)),
         S_T3 = mean(c(S_T3_rep1, S_T3_rep2, S_T3_rep3)),
         S_T4 = mean(c(S_T4_rep1, S_T4_rep2, S_T4_rep3)),
         .keep = "unused")

str(sctc_cpm)
# rowwise_df [1,293 × 9] (S3: rowwise_df/tbl_df/tbl/data.frame)

write.table(
  sctc_cpm,
  "results/SCTC.rowmean_cpm.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
