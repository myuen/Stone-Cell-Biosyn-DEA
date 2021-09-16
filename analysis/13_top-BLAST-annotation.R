library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)


### Read BLAST results for all DE contigs
blast <- read.table("data/SCB.upRegDSC.17Jun.blastpTAIR10.txt", 
                    header = TRUE, sep = "\t", 
                    quote = "", stringsAsFactors = FALSE)

blast.tophit <- blast %>% nest(data = -qseqid)

blast.tophit$data <- blast.tophit$data %>%
  map(function(x){
  x[1,]})

blast.tophit <- blast.tophit %>% unnest(-qseqid)

blast.tophit <- blast.tophit %>% select(qseqid, sseqid, salltitles)

blast.tophit$salltitles <- str_replace(blast.tophit$salltitles, "\\W+Symbols:\\s", "")

str(blast.tophit)
# tibble [1,025 × 3] (S3: tbl_df/tbl/data.frame)


# Genes of interest
goi <- c("NST", "VND", "VNS", "MYB")

goi <- 
  map_dfr(goi, function(g){
  blast.tophit[str_detect(blast.tophit$salltitles, g),]
})

str(goi)
# A tibble: 23 x 3


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb.cpm <- read.delim("results/SCB.tmm_normalized_cpm.txt")

scb.cpm <- scb.cpm %>%
  rownames_to_column("cds") %>%
  filter(cds %in% blast.tophit$qseqid)

str(scb.cpm)
# 'data.frame':	1025 obs. of  13 variables:

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
# rowwise_df [1,025 × 5] (S3: rowwise_df/tbl_df/tbl/data.frame)


# Read TMM normalized CPM mapped to differentially expressed CDS in 
# LMD experiment
sctc.cpm <- read.delim("results/SCTC.tmm_normalized_cpm.txt")

sctc.cpm <- sctc.cpm %>%
  rownames_to_column("cds") %>%
  filter(cds %in% blast.tophit$qseqid)

str(sctc.cpm)
# 'data.frame':	1025 obs. of  25 variables:

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

str(sctc.cpm)
# rowwise_df [1,025 × 9] (S3: rowwise_df/tbl_df/tbl/data.frame)

###

scb.goi <- inner_join(goi, scb.cpm, by = c("qseqid" = "cds"))

scb.goi <- scb.goi %>% 
  column_to_rownames("qseqid")

str(scb.goi)
# 'data.frame':	23 obs. of  6 variables:

###

sctc.goi <- inner_join(goi, sctc.cpm, by = c("qseqid" = "cds"))

sctc.goi <- sctc.goi %>% 
  column_to_rownames("qseqid")

str(sctc.goi)
# 'data.frame':	23 obs. of  10 variables:

sctc.r.goi <- sctc.goi[,!c(str_detect(colnames(sctc.goi), "^S"))]
sctc.r.goi <- sctc.r.goi[rowSums(sctc.r.goi[,3:6]) != 0,]

sctc.s.goi <- sctc.goi[,!c(str_detect(colnames(sctc.goi), "^R"))]
sctc.s.goi <- sctc.s.goi[rowSums(sctc.s.goi[,3:6]) != 0,]


makeHM <- function(goi, out_prefix) {

  goi.cpm <- goi %>% 
    select(-c(sseqid, salltitles)) %>% as.matrix()
  hm.dist <- as.dist(1- cor(t(goi.cpm), method = "pearson"))
  hm.clust <- hclust(hm.dist, method = "complete")
  hm.dendrogram <- as.dendrogram(hm.clust)
  
  annot <- goi$salltitles %>% 
    str_replace("\\W+Symbols:\\s", "")

  hovertext <- goi %>% select(-c(sseqid, salltitles))
  
  hovertext[, 1:dim(hovertext)[2]] <- annot
  
  hm <- heatmaply(goi.cpm,
                  colors = heat.colors(50),
                  scale = "row",
                  
                  Rowv = hm.dendrogram,
                  show_dendrogram = c(TRUE, FALSE),
                  row_dend_left = TRUE,
                  
                  hide_colorbar = TRUE,
                  dendrogram = "row",
                  showticklabels = c(TRUE, FALSE),
                  
                  custom_hovertext = hovertext
  )
  
  htmlwidgets::saveWidget(
    hm,
    file = paste("results/figures/", out_prefix, ".html", sep = ""),
    selfcontained = TRUE)

  # return (hm)
}

makeHM(scb.goi, "lmd-mockup-hm.sep16")
makeHM(sctc.r.goi, "sctc-r-mockup-hm.sep16")
makeHM(sctc.s.goi, "sctc-s-mockup-hm.sep16")


### END HERE ###


sctc.s_peakedT2 <- 
  scan("results/SCTC.S_PeakedT2.txt", what = "character")
# Read 31 items

sctc.s_peakedT3 <- 
  scan("results/SCTC.S_PeakedT3.txt", what = "character")
# Read 26 items

sctc.r_peakedT2 <- 
  scan("results/SCTC.R_PeakedT2.txt", what = "character")
# Read 158 items

sctc.r_peakedT3 <- 
  scan("results/SCTC.R_PeakedT3.txt", what = "character")
# Read 342 items


### S T2
sctc.s_peakedT2.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.s_peakedT2)

str(sctc.s_peakedT2.tophit)
# tibble [18 × 21] (S3: tbl_df/tbl/data.frame)

# write.table(sctc.s_peakedT2.tophit, 
#             "results/SCTC.S_PeakedT2.blastTopHit.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

write.table(sctc.s_peakedT2.tophit, 
            "results/SCTC.S_PeakedT2.blastTopTAIR10Hit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### S T3
sctc.s_peakedT3.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.s_peakedT3)

str(sctc.s_peakedT3.tophit)
# tibble [20 × 21] (S3: tbl_df/tbl/data.frame)

# write.table(sctc.s_peakedT3.tophit, 
#             "results/SCTC.S_PeakedT3.blastTopHit.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

write.table(sctc.s_peakedT3.tophit,
            "results/SCTC.S_PeakedT3.blastTopTAIR10Hit.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T2
sctc.r_peakedT2.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.r_peakedT2)

str(sctc.r_peakedT2.tophit)
# tibble [121 × 21] (S3: tbl_df/tbl/data.frame)

# str(sctc.r_peakedT2.tophit)
# tibble [143 × 21] (S3: tbl_df/tbl/data.frame)

# write.table(sctc.r_peakedT2.tophit, 
#             "results/SCTC.R_PeakedT2.blastTopHit.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

write.table(sctc.r_peakedT2.tophit, 
            "results/SCTC.R_PeakedT2.blastTopTAIR10Hit.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T3
sctc.r_peakedT3.tophit <-
  blast.tophit %>%
  filter(qseqid %in% sctc.r_peakedT3)

str(sctc.r_peakedT3.tophit)
# tibble [261 × 21] (S3: tbl_df/tbl/data.frame)

# write.table(sctc.r_peakedT3.tophit, 
#             "results/SCTC.R_PeakedT3.blastTopHit.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

write.table(sctc.r_peakedT3.tophit,
            "results/SCTC.R_PeakedT3.blastTopTAIR10Hit.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")
