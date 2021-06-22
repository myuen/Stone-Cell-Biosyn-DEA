library(dplyr)
library(heatmaply)
# library(tidyr)
library(stringr)
library(tibble)
library(tidyselect)


# Read statistics from differential expression analysis for LMD experiment
sigDE <- read.delim("results/SCB.sigDE_stats.17Jun.txt",
                    stringsAsFactors = FALSE,
                    header = TRUE)
str(sigDE)
# 'data.frame':	2283 obs. of  4 variables:


# We are only interested in those that are up-regulated in 
# developing stone cells (DSC)
upReg <- sigDE %>% filter(logFC >= 0)


# Identify those that are up-regulated only in Q903
S_cType_excl <- upReg %>% 
  filter(focus_term == "S_cType" & focus_term != "R_cType") %>%
  select("cds") %>%
  distinct()

str(S_cType_excl)
# 'data.frame':	244 obs. of  1 variable:

# Identify those that are up-regulated only in H898
R_cType_excl <- upReg %>% 
  filter(focus_term != "S_cType" & focus_term == "R_cType") %>%
  select("cds") %>%
  distinct()

str(R_cType_excl)
# 'data.frame':	1129 obs. of  1 variable:


# Identify those that are up-regulated in both genotypes 
commons <- upReg[duplicated(upReg$cds), "cds"]

commons <- unique(commons)

str(commons)
# chr [1:80]


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb <- read.delim("results/SCB.tmm_normalized_cpm.txt")

scb <- scb %>%
  rownames_to_column("cds")


# Can't get mutate or transmute to work with select(starts_with())
scb <- scb %>%
  rowwise() %>%
  mutate(R_CP = mean(c(R_CP_rep1, R_CP_rep2, R_CP_rep3)))

scb <- scb %>%
  rowwise() %>%
  mutate(R_DSC = mean(c(R_DSC_rep1, R_DSC_rep2, R_DSC_rep3)))

scb <- scb %>%
  rowwise() %>%
  mutate(S_CP = mean(c(S_CP_rep1, S_CP_rep2, S_CP_rep3)))

scb <- scb %>%
  rowwise() %>%
  mutate(S_DSC = mean(c(S_DSC_rep1, S_DSC_rep2, S_DSC_rep3)))

scb <- scb %>% 
  select(-matches("rep\\d"))

str(scb)
# rowwise_df [26,787 × 5] (S3: rowwise_df/tbl_df/tbl/data.frame)



# Read TMM normalized CPM from mapped to differentially expressed CDS in 
# LMD experiment
sctc <- read.delim("results/SCTC.tmm_normalized_cpm.txt")

sctc <- sctc %>% rownames_to_column("cds")

str(sctc)
# 'data.frame':	5305 obs. of  25 variables:

sctc <- sctc %>%
  rowwise() %>%
  mutate(R_T1 = mean(c(R_T1_rep1, R_T1_rep2, R_T1_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(R_T2 = mean(c(R_T2_rep1, R_T2_rep2, R_T2_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(R_T3 = mean(c(R_T3_rep1, R_T3_rep2, R_T3_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(R_T4 = mean(c(R_T4_rep1, R_T4_rep2, R_T4_rep3)))


sctc <- sctc %>%
  rowwise() %>%
  mutate(S_T1 = mean(c(S_T1_rep1, S_T1_rep2, S_T1_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(S_T2 = mean(c(S_T2_rep1, S_T2_rep2, S_T2_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(S_T3 = mean(c(S_T3_rep1, S_T3_rep2, S_T3_rep3)))

sctc <- sctc %>%
  rowwise() %>%
  mutate(S_T4 = mean(c(S_T4_rep1, S_T4_rep2, S_T4_rep3)))

sctc <- sctc %>% 
  select(-matches("_rep\\d"))

str(sctc)
# rowwise_df [5,305 × 9] (S3: rowwise_df/tbl_df/tbl/data.frame)


# Joining two dataframes.  These are CPM for all DE contigs in LMD experiment
sigDE_cpm <- right_join(scb, sctc)
# Joining, by = "cds"

str(sigDE_cpm)
# rowwise_df [5,305 × 13] (S3: rowwise_df/tbl_df/tbl/data.frame)


# Subset sigDE_cpm for S_cType_excl
S_cType_excl_cpm <- sigDE_cpm %>% 
  filter(cds %in% S_cType_excl$cds)

str(S_cType_excl_cpm)
# rowwise_df [244 × 13] (S3: rowwise_df/tbl_df/tbl/data.frame)


# Subset sigDE_cpm for R_cType_excl
R_cType_excl_cpm <- sigDE_cpm %>%
  filter(cds %in% R_cType_excl$cds)

str(R_cType_excl_cpm)
# rowwise_df [1,129 × 13] (S3: rowwise_df/tbl_df/tbl/data.frame)




# Takes 2 matrices and 2 output file as argument
# Heatmap 2 is base on distrance matrix created from heatmap 1
makeHM <- function(scb_cpm, sctc_cpm, outfile) {
  
  scb_cpm <- as.matrix(scb_cpm)
  
  scb_d <- as.dist(1- cor(t(scb_cpm), method = "pearson"))
  
  scb_c <- hclust(scb_d, method = "complete")

  scb_dendrogram <- as.dendrogram(scb_c)
  
  scb_hm <- 
    heatmaply(scb_cpm,
              colors = heat.colors(50),
              scale = "row",
              row_dend_left = TRUE,
              Rowv = scb_dendrogram,
              hide_colorbar = TRUE,
              dendrogram = "row",
              showticklabels = c(TRUE, FALSE)
              # file = outfile1
              )
  
  sctc_hm <- 
    heatmaply(as.matrix(sctc_cpm),
              colors = heat.colors(50),
              scale = "row",
              Rowv = scb_dendrogram,
              hide_colorbar = TRUE,
              dendrogram = "row",
              show_dendrogram = c(FALSE, FALSE),
              showticklabels = c(TRUE, FALSE)
              # file = outfile2
              )
 
  # merged <- subplot(scb_hm, sctc_hm, margin = 0.005)
  merged <- subplot(scb_hm, sctc_hm, margin = 0)
  
  htmlwidgets::saveWidget(merged, 
                          file = outfile, 
                          selfcontained = TRUE)
}

# Subset sigDE_cpm for commons
commons_cpm <- sigDE_cpm[sigDE_cpm$cds %in% commons, ]

str(commons_cpm)
# rowwise_df [80 × 13] (S3: rowwise_df/tbl_df/tbl/data.frame)

scb.commons <- scb %>% 
  filter(cds %in% commons) %>%
  column_to_rownames("cds")

sctc.commons <- sctc %>%
  filter(cds %in% commons) %>%
  column_to_rownames("cds")

makeHM(scb.commons, sctc.commons, "results/figures/commons.html")


# Subset sigDE_cpm for S_cType_excl
scb.s_ctype <- scb %>% 
  filter(cds %in% S_cType_excl$cds) %>%
  column_to_rownames("cds")

sctc.s_ctype <- sctc %>% 
  filter(cds %in% S_cType_excl$cds) %>%
  column_to_rownames("cds")

makeHM(scb.s_ctype, sctc.s_ctype, "results/figures/S_cType.html")


# Subset sigDE_cpm for R_cType_excl
scb.r_ctype <- scb %>% 
  filter(cds %in% R_cType_excl$cds) %>%
  column_to_rownames("cds")

sctc.r_ctype <- sctc %>% 
  filter(cds %in% R_cType_excl$cds) %>%
  column_to_rownames("cds")

makeHM(scb.r_ctype, sctc.r_ctype, "results/figures/R_cType.html")




# makeHM(S_cType_excl_cpm, "results/figures/S_cType_hm.html")
# 
# makeHM(R_cType_excl_cpm, "results/figures/R_cType_hm.html")
# 
# makeHM(commons_cpm, "results/figures/commons.html")
