library(dplyr)
library(heatmaply)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

library(org.At.tair.db)
columns(org.At.tair.db)


secondary_cellwall_goIDs <- 
  c("GO:0030244", # Cellulose biosynthesis genes
    "GO:0009834", # Secondary cel wall biosynthesis
    "GO:0009809"  # Lignin biosynthesis genes
  )

go_2_tair <- as.list(org.At.tairGO2ALLTAIRS)

go_2_tair <- go_2_tair[secondary_cellwall_goIDs]

# tairids <- map(goid_tairids, function(g){
#   # g[names(g) == "TAS"]
#   g[names(g) %in% c("IC","TAS")]
# }) %>% unlist() %>% unique()
# 
# length(tairids)
# [1] 11


# TAIR IDs associated with secondary cell wall biosynthesis
secondary_cellwall_tairIDs <- go_2_tair %>% unlist() %>% unique()

length(secondary_cellwall_tairIDs)
# [1] 183


### Read BLAST results for all DE contigs
blast <- read.table("data/SCB.upRegDSC.17Jun.blastpTAIR10.txt", 
                    header = TRUE, sep = "\t", 
                    quote = "", stringsAsFactors = FALSE)

blast.tophit <- blast %>% nest(data = -qseqid)

blast.tophit$data <- blast.tophit$data %>%
  map(function(x){
    x[1,]})

blast.tophit <- blast.tophit %>% unnest(-qseqid)

blast.tophit <- blast.tophit %>% 
  dplyr::select(qseqid, sseqid, salltitles)

blast.tophit$salltitles <- 
  str_replace(blast.tophit$salltitles, "\\W+Symbols:\\s", "")

# Strip TAIR10 version ID from accession number
blast.tophit$sseqid <- 
  blast.tophit$sseqid %>% 
  str_replace(".\\d$", "")

# BLAST hit pertaining to secondary cell wall biosynthesis
blast.tophit <- blast.tophit %>% 
  filter(sseqid %in% secondary_cellwall_tairIDs)

str(blast.tophit)
# tibble [73 × 3] (S3: tbl_df/tbl/data.frame)



# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
scb_rowmean_cpm <- read.delim("results/SCB.rowmean_cpm.txt")

scb.2ndCW <- 
  inner_join(blast.tophit, scb_rowmean_cpm, by = c("qseqid" = "cds")) %>%
  column_to_rownames("qseqid")

str(scb.2ndCW)
# 'data.frame':	73 obs. of  6 variables:


# Read TMM normalized CPM mapped to differentially expressed CDS in 
# LMD experiment
sctc_rowmean_cpm <- read.delim("results/SCTC.rowmean_cpm.txt")

sctc.2ndCW <-
  inner_join(blast.tophit, sctc_rowmean_cpm, by = c("qseqid" = "cds")) %>%
  column_to_rownames("qseqid")

str(sctc.2ndCW)
# 'data.frame':	73 obs. of  10 variables:



sctc.r.2ndCW <- sctc.2ndCW %>% 
  dplyr::select(!starts_with("S", ignore.case = FALSE))

# TODO: Check CDS expression != 0
# sctc.r.2ndCW <- sctc.r.2ndCW[,rowSums(c(str_starts(colnames(sctc.r.2ndCW), "R")))]

str(sctc.r.2ndCW)
# 'data.frame':	73 obs. of  6 variables:


sctc.s.2ndCW <- sctc.2ndCW %>% 
  dplyr::select(!starts_with("R", ignore.case = FALSE))

# TODO: Check CDS expression != 0
# sctc.s.2ndCW <- 
#   sctc.r.2ndCW %>% dplyr::select(starts_with("R")) %>% rowSums()

str(sctc.s.2ndCW)
# 'data.frame':	73 obs. of  6 variables:

###

makeHM <- function(goi, out_prefix) {
  
  goi.cpm <- goi %>% 
    dplyr::select(-c(sseqid, salltitles)) %>% as.matrix()
  hm.dist <- as.dist(1- cor(t(goi.cpm), method = "pearson"))
  hm.clust <- hclust(hm.dist, method = "complete")
  hm.dendrogram <- as.dendrogram(hm.clust)
  
  annot <- goi$salltitles %>% 
    str_replace("\\W+Symbols:\\s", "")
  
  hovertext <- goi %>% dplyr::select(-c(sseqid, salltitles))
  
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

makeHM(scb.2ndCW, "lmd-2nd-cellwall-hm.sep23")
makeHM(sctc.r.2ndCW, "sctc-r-2nd-cellwall-hm.sep23")
makeHM(sctc.s.2ndCW, "sctc-s-2nd-cellwall-hm.sep23")
