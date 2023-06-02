library(edgeR)
library(limma)
library(readr)
library(stringr)
library(tximport)


source("analysis/helper02_TC-PCA-maker.R")


# Import file with tximport
quant.files <-
  dir("data/SCTC_Salmon-17Jun",
      pattern = "quant.sf",
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE)

samples <-
  str_extract(quant.files, "\\d{3}_\\w+")
  
samples <- samples %>%
  str_replace("-", ".")

samples <- samples %>%
  str_replace("898", "R") %>%
  str_replace("903", "S")

samples <- samples %>%
  str_replace("A6", "T1") %>%
  str_replace("A7", "T2") %>%
  str_replace("A8", "T3") %>%
  str_replace("A10", "T4")

samples <- paste(samples, "_rep", sep = "") %>% 
  paste(rep(c(1,2,3), 8), sep = "")

names(quant.files) <- samples


txi <-
  tximport(quant.files,
           type = "salmon",
           txOut = TRUE,
           varReduce = TRUE,
           countsFromAbundance = "lengthScaledTPM")


# Create experimental design
expDes <- 
  data.frame(
    sample = samples,
    gType = factor(c(samples %>% str_extract("^\\w")),
                   levels = c("S", "R")),
    cType = factor(c(samples %>% 
                       str_extract("_\\w\\d_") %>% 
                       str_replace_all("_", "")),
                   levels = c("T1", "T2", "T3", "T4")),
    bioRep = as.numeric(samples %>% str_extract("\\d$")))

expDes$group <- with(expDes, interaction(gType, cType))

#       sample gType cType bioRep group
# 1  R_T4_rep1     R    T4      1  R.T4
# 2  R_T4_rep2     R    T4      2  R.T4
# 3  R_T4_rep3     R    T4      3  R.T4
# 4  R_T1_rep1     R    T1      1  R.T1
# 5  R_T1_rep2     R    T1      2  R.T1
# 6  R_T1_rep3     R    T1      3  R.T1
# 7  R_T2_rep1     R    T2      1  R.T2
# 8  R_T2_rep2     R    T2      2  R.T2
# 9  R_T2_rep3     R    T2      3  R.T2
# 10 R_T3_rep1     R    T3      1  R.T3
# 11 R_T3_rep2     R    T3      2  R.T3
# 12 R_T3_rep3     R    T3      3  R.T3
# 13 S_T4_rep1     S    T4      1  S.T4
# 14 S_T4_rep2     S    T4      2  S.T4
# 15 S_T4_rep3     S    T4      3  S.T4
# 16 S_T1_rep1     S    T1      1  S.T1
# 17 S_T1_rep2     S    T1      2  S.T1
# 18 S_T1_rep3     S    T1      3  S.T1
# 19 S_T2_rep1     S    T2      1  S.T2
# 20 S_T2_rep2     S    T2      2  S.T2
# 21 S_T2_rep3     S    T2      3  S.T2
# 22 S_T3_rep1     S    T3      1  S.T3
# 23 S_T3_rep2     S    T3      2  S.T3
# 24 S_T3_rep3     S    T3      3  S.T3


# Load counts into DGEList object from edgeR package.
x <- DGEList(counts = txi$counts, 
             group = expDes$group)

dim(x)
# [1] 1293   24


# TMM Normalization by Depth
y <- calcNormFactors(x)

# tc for TimeCourse
tc.cpm <- cpm(y)

write.table(
  sctc.cpm,
  "results/TC.tmm_normalized_cpm.txt",
  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

v <- voom(y, plot = TRUE)


(p <- PCA_maker(expDes, v))

ggsave("results/figures/TC-PCA.svg", plot = p,
       height = 6, width = 6)
