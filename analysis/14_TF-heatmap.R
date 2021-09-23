library(dplyr)
library(heatmaply)
library(purrr)
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

blast.tophit$salltitles <- 
  str_replace(blast.tophit$salltitles, "\\W+Symbols:\\s", "")

str(blast.tophit)
# tibble [1,025 × 3] (S3: tbl_df/tbl/data.frame)


# transcription factors of interest
tf <- c("NST", "VND", "VNS", "MYB")

tf <- 
  map_dfr(tf, function(f){
  blast.tophit[str_detect(blast.tophit$salltitles, f),]
})

str(tf)
# A tibble: 23 x 3


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb_rowmean_cpm <- read.delim("results/SCB.rowmean_cpm.txt")

scb_rowmean_cpm <- scb_rowmean_cpm %>%
  filter(cds %in% blast.tophit$qseqid)

str(scb_rowmean_cpm)
# 'data.frame':	1025 obs. of  5 variables:


# Read TMM normalized CPM mapped to differentially expressed CDS in 
# LMD experiment
sctc_rowmean_cpm <- read.delim("results/SCTC.rowmean_cpm.txt")

sctc_rowmean_cpm <- sctc_rowmean_cpm %>%
  filter(cds %in% blast.tophit$qseqid)

str(sctc_rowmean_cpm)
# 'data.frame':	1025 obs. of  25 variables:


scb.tf <- inner_join(tf, scb_rowmean_cpm, by = c("qseqid" = "cds"))

scb.tf <- scb.tf %>% 
  column_to_rownames("qseqid")

str(scb.tf)
# 'data.frame':	23 obs. of  6 variables:

###

sctc.tf <- inner_join(tf, sctc_rowmean_cpm, by = c("qseqid" = "cds"))

sctc.tf <- sctc.tf %>% 
  column_to_rownames("qseqid")

str(sctc.tf)
# 'data.frame':	23 obs. of  10 variables:

sctc.r.tf <- sctc.tf[,!c(str_detect(colnames(sctc.tf), "^S"))]
sctc.r.tf <- sctc.r.tf[rowSums(sctc.r.tf[,3:6]) != 0,]
str(sctc.r.tf)
# 'data.frame':	23 obs. of  6 variables:

sctc.s.tf <- sctc.tf[,!c(str_detect(colnames(sctc.tf), "^R"))]
sctc.s.tf <- sctc.s.tf[rowSums(sctc.s.tf[,3:6]) != 0,]
str(sctc.s.tf)
# 'data.frame':	22 obs. of  6 variables:


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

makeHM(scb.tf, "lmd-tf-hm.sep23")
makeHM(sctc.r.tf, "sctc-r-tf-hm.sep23")
makeHM(sctc.s.tf, "sctc-s-tf-hm.sep23")
