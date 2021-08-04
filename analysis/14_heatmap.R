library(dplyr)
library(heatmaply)
# library(tidyr)
# library(orca)
library(stringr)
library(tibble)
library(tidyselect)

source("analysis/helper02_heatmap-maker.R")


# Read statistics from differential expression analysis for LMD experiment
sigDE <- read.delim("results/SCB.sigDE_stats.17Jun.txt",
                    stringsAsFactors = FALSE,
                    header = TRUE)

str(sigDE)
# 'data.frame':	2283 obs. of  4 variables:


# We are only interested in those that are up-regulated in 
# developing stone cells (DSC)
scb_upReg <- sigDE %>% filter(logFC >= 0)

str(scb_upReg)
# 'data.frame':	1373 obs. of  4 variables:


# Identify those that are up-regulated only in susceptible genotype
scb.s_upReg <- scb_upReg %>% 
  filter(focus_term == "S_cType" & focus_term != "R_cType") %>%
  select("cds") %>%
  distinct()

scb.s_upReg <- as.character(scb.s_upReg$cds)

str(scb.s_upReg)
#  chr [1:244] 


# Identify those that are up-regulated only in resistance genotype
scb.r_upReg <- scb_upReg %>% 
  filter(focus_term != "S_cType" & focus_term == "R_cType") %>%
  select("cds") %>%
  distinct()

scb.r_upReg <- as.character(scb.r_upReg$cds)

str(scb.r_upReg)
#  chr [1:1129]


# Identify those that are up-regulated in both genotypes
scb.rs_upReg <- scb_upReg[duplicated(scb_upReg$cds), "cds"]

str(scb.rs_upReg)
# chr [1:80]


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb.cpm <- read.delim("results/SCB.tmm_normalized_cpm.txt")

scb.cpm <- scb.cpm %>%
  rownames_to_column("cds") %>%
  filter(cds %in% scb_upReg$cds)


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
  filter(cds %in% scb_upReg$cds)

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

scb.s_upReg.hm <- makeHM(scb.s_upReg.cpm, "something")

orca(scb.s_upReg.hm,
     file = "results/figures/scb.s_upReg.9Jul.svg")


# 1b - 
sctc.s_peakedT2T3.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.s_upReg) %>%
  filter((S_T2 > S_T1) & (S_T2 > S_T4) & (S_T3 > S_T1) & (S_T3 > S_T4)) %>%
  column_to_rownames("cds")

str(sctc.s_peakedT2T3.cpm)
# 'data.frame':	34 obs. of  8 variables:

write(rownames(sctc.s_peakedT2T3.cpm), 
          "data/sctc.s_peakedT2T3.cdsId.txt")

sctc.s_peakedT2T3.hm <- 
  make2HM(sctc.s_peakedT2T3.cpm %>% select(starts_with("S_T")),
          sctc.s_peakedT2T3.cpm %>% select(starts_with("R_T")),
          "sctc.s_peakedT2T3")

orca(sctc.s.peakedT2T3.hm,
     file = "results/figures/sctc.s_peakedT2T3.9Jul.svg")


# 2a - Subset LMD data for sequences up-regulated in resistance genotype only
scb.r_upReg.cpm <- scb.cpm %>%
  filter(cds %in% scb.r_upReg) %>%
  column_to_rownames("cds")

str(scb.r_upReg.cpm)
# 'data.frame':	1129 obs. of  4 variables:

scb.r_upReg.hm <- makeHM(scb.r_upReg.cpm, "something")

orca(scb.r_upReg.hm,
     file = "results/figures/scb.r_upReg.9Jul.svg")


# 2b
sctc.r_peakedT2T3.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.r_upReg) %>%
    filter((R_T2 > R_T1) & (R_T2 > R_T4) & (R_T3 >= R_T1) & (R_T3 > R_T4)) %>%
  column_to_rownames("cds")

str(sctc.r_peakedT2T3.cpm)
# 'data.frame':	256 obs. of  8 variables:

write(rownames(sctc.r_peakedT2T3.cpm), 
      "data/sctc.r_peakedT2T3.cdsId.txt")

sctc.r_peakedT2T3.hm <- 
  make2HM(sctc.r_peakedT2T3.cpm %>% select(starts_with("R_T")),
          sctc.r_peakedT2T3.cpm %>% select(starts_with("S_T")),
        "sctc.r_peakedT2T3")

orca(sctc.r_peakedT2T3.hm,
     file = "results/figures/sctc.r_peakedT2T3.9Jul.svg")


# 3 - Subset sigDE_cpm for contig up-regulated in both resistance and 
# susceptible genotype
scb.rs_upReg.cpm <- scb.cpm %>%
  filter(cds %in% scb.rs_upReg) %>%
  column_to_rownames("cds")

str(scb.rs_upReg.cpm)
# 'data.frame':	80 obs. of  4 variables:

scb.rs_upReg.hm <- makeHM(scb.rs_upReg.cpm, "something")

orca(scb.rs_upReg.hm,
     file = "results/figures/scb.rs_upReg.9Jul.svg")

# 3b
sctc.rs_peakedT2T3.cpm <- sctc.cpm %>% 
  filter(cds %in% scb.rs.upReg) %>%
  filter((R_T2 > R_T1) & (R_T2 > R_T4) & (R_T3 >= R_T1) & (R_T3 > R_T4)) %>%
  column_to_rownames("cds")

str(sctc.rs_peakedT2T3.cpm)
# 'data.frame':	19 obs. of  8 variables:

sctc.rs_peakedT2T3.hm <-   
  make2HM(sctc.rs_peakedT2T3.cpm %>% select(starts_with("R_T")),
          sctc.rs_peakedT2T3.cpm %>% select(starts_with("S_T")),
          "sctc.rs_peakedT2T3")

orca(sctc.rs_peakedT2T3.hm,
     file = "results/figures/sctc.rs_peakedT2T3.9Jul.svg")
