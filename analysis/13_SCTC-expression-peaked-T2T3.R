library(dplyr)
library(tibble)


# Read statistics from differential expression analysis for LMD experiment
sigDE <- read.delim("results/SCB.sigDE_stats.17Jun.txt",
                    stringsAsFactors = FALSE,
                    header = TRUE)

str(sigDE)
# 'data.frame':	2283 obs. of  4 variables:


# We are only interested in those that are up-regulated in 
# developing stone cells (DSC)
scb.upReg <- sigDE %>% filter(logFC >= 0)

str(scb.upReg)
# 'data.frame':	1373 obs. of  4 variables:


# Identify those that are up-regulated only in susceptible genotype
scb.s_upReg <- scb.upReg %>% 
  filter(focus_term == "S_cType" & focus_term != "R_cType") %>%
  select("cds") %>%
  distinct()

scb.s_upReg <- as.character(scb.s_upReg$cds)

str(scb.s_upReg)
#  chr [1:244] 


# Identify those that are up-regulated only in resistance genotype
scb.r_upReg <- scb.upReg %>% 
  filter(focus_term != "S_cType" & focus_term == "R_cType") %>%
  select("cds") %>%
  distinct()

scb.r_upReg <- as.character(scb.r_upReg$cds)

str(scb.r_upReg)
#  chr [1:1129]


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb.cpm <- read.delim("results/SCB.tmm_normalized_cpm.txt")

scb.cpm <- scb.cpm %>%
  rownames_to_column("cds") %>%
  filter(cds %in% scb.upReg$cds)


# This is cumbersome but cannot get mutate to work 
# with rowMeans(select(starts_with(col)))
scb.cpm <- scb.cpm %>% 
  rowwise() %>% 
  mutate(R_CP = mean(c(R_CP_rep1, R_CP_rep2, R_CP_rep3)),
         R_DSC = mean(c(R_DSC_rep1, R_DSC_rep2, R_DSC_rep3)),
         S_CP = mean(c(S_CP_rep1, S_CP_rep2, S_CP_rep3)),
         S_DSC = mean(c(S_DSC_rep1, S_DSC_rep2, S_DSC_rep3)),
         .keep = "unused")

str(scb.cpm)
# rowwise_df [1,293 × 5] (S3: rowwise_df/tbl_df/tbl/data.frame)


# Read TMM normalized CPM mapped to differentially expressed CDS in 
# LMD experiment
sctc.cpm <- read.delim("results/SCTC.tmm_normalized_cpm.txt")

sctc.cpm <- sctc.cpm %>% 
  rownames_to_column("cds") %>%
  filter(cds %in% scb.upReg$cds)

sctc.cpm <- sctc.cpm %>%
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

sctc.cpm <- sctc.cpm %>% 
  select(-matches("_rep\\d"))

str(sctc.cpm)
# rowwise_df [1,293 × 8] (S3: rowwise_df/tbl_df/tbl/data.frame)


# 1a - Subset LMD data for sequences up-regulated in susceptible genotype only
scb.s_upReg.cpm <- scb.cpm %>%
  filter(cds %in% scb.s_upReg) %>%
  column_to_rownames("cds")

str(scb.s_upReg.cpm)
# 'data.frame':	244 obs. of  4 variables:


# 1b - Subset temporal data for sequences 
# up-regulated and expression peaked in T2
sctc.s_peakedT2.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.s_upReg) %>%
  filter((S_T2 > S_T1) & (S_T2 > S_T4) & (S_T2 > S_T3)) %>%
  column_to_rownames("cds")

str(sctc.s_peakedT2.cpm)
# 'data.frame':	31 obs. of  8 variables:

write(rownames(sctc.s_peakedT2.cpm), 
            "results/SCTC.S_PeakedT2.txt")


# 1c - Subset temporal data for sequences
# up-regulated and expression peaked in T3
sctc.s_peakedT3.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.s_upReg) %>%
  filter((S_T3 > S_T1) & (S_T3 > S_T4) & (S_T3 > S_T2)) %>%
  column_to_rownames("cds")

str(sctc.s_peakedT3.cpm)
# 'data.frame':	26 obs. of  8 variables:

write(rownames(sctc.s_peakedT3.cpm),
      "results/SCTC.S_PeakedT3.txt")


# 2a - Subset LMD data for sequences up-regulated in resistance genotype only
scb.r_upReg.cpm <- scb.cpm %>%
  filter(cds %in% scb.r_upReg) %>%
  column_to_rownames("cds")

str(scb.r_upReg.cpm)
# 'data.frame':	1129 obs. of  4 variables:


# 2b - Subset temporal data for sequences 
# up-regulated and expression peaked in T2
sctc.r_peakedT2.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.r_upReg) %>%
  filter((R_T2 > R_T1) & (R_T2 > R_T4) & (R_T2 > R_T3)) %>%
  column_to_rownames("cds")

str(sctc.r_peakedT2.cpm)
# 'data.frame':	158 obs. of  8 variables:

write(rownames(sctc.r_peakedT2.cpm), 
      "results/SCTC.R_PeakedT2.txt")


# 2c - Subset temporal data for sequences 
# up-regulated and expression peaked in T3
sctc.r_peakedT3.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.r_upReg) %>%
  filter((R_T3 > R_T1) & (R_T3 > R_T4) & (R_T3 > R_T2)) %>%
  column_to_rownames("cds")

str(sctc.r_peakedT3.cpm)
# 'data.frame':	342 obs. of  8 variables:

write(rownames(sctc.r_peakedT3.cpm), 
      "results/SCTC.R_PeakedT3.txt")
