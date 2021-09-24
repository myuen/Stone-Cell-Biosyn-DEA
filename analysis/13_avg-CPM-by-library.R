library(dplyr)
library(tibble)

options(digits = 15)

# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
lmd.cpm <- read.delim("results/LMD.tmm_normalized_cpm.txt")

lmd.cpm <- lmd.cpm %>% 
  transmute("R_CP" = rowMeans(select(lmd.cpm, starts_with("R_CP"))),
            "R_DSC" = rowMeans(select(lmd.cpm, starts_with("R_DSC"))),
            "S_CP" = rowMeans(select(lmd.cpm, starts_with("S_CP"))),
            "S_DSC" = rowMeans(select(lmd.cpm, starts_with("S_DSC")))) %>%
  rownames_to_column("cds")

str(lmd.cpm)
# 'data.frame':	26787 obs. of  5 variables:

write.table(
  lmd.cpm,
  "results/LMD.rowmean_cpm.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# Read TMM normalized CPM from time-course (TC) experiment mapped to 
# differentially expressed CDS in LMD experiment
tc.cpm <- read.delim("results/TC.tmm_normalized_cpm.txt")

tc.cpm <- tc.cpm %>%
  transmute(R_T1 = rowMeans(select(tc.cpm, starts_with("R_T1"))),
            R_T2 = rowMeans(select(tc.cpm, starts_with("R_T2"))),
            R_T3 = rowMeans(select(tc.cpm, starts_with("R_T3"))),
            R_T4 = rowMeans(select(tc.cpm, starts_with("R_T4"))),
            S_T1 = rowMeans(select(tc.cpm, starts_with("S_T1"))),
            S_T2 = rowMeans(select(tc.cpm, starts_with("S_T2"))),
            S_T3 = rowMeans(select(tc.cpm, starts_with("S_T3"))),
            S_T4 = rowMeans(select(tc.cpm, starts_with("S_T4")))) %>% 
  rownames_to_column("cds")

str(tc.cpm)
# 'data.frame':	1293 obs. of  9 variables:

write.table(
  tc.cpm,
  "results/TC.rowmean_cpm.txt",
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
