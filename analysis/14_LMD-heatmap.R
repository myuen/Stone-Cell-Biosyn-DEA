library(dplyr)
library(heatmaply)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

source("analysis/helper02_heatmap-maker.R")

# Read CDS ID that are up-reguled in the LMD experiment
lmd.upReg <- scan("results/LMD.upRegDSC.txt", what = "character")
# Read 1293 items


# Read TMM normalized CPM from differentially expressed CDS from LMD experiment
lmd.rowmean.cpm <- read.delim("results/LMD.rowmean_cpm.txt")

str(lmd.rowmean.cpm)
# 'data.frame':	26787 obs. of  5 variables:

lmd.upReg.cpm <- lmd.rowmean.cpm %>% 
  filter(cds %in% lmd.upReg) %>% 
  column_to_rownames("cds")

str(lmd.upReg.cpm)
# 'data.frame':	1293 obs. of  5 variables:


makeHM(lmd.upReg.cpm, "lmd-hm")

