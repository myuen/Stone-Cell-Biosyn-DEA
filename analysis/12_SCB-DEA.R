require(edgeR)
require(dplyr)
require(limma)
require(purrr)
require(tibble)
require(tximport)

# source("analysis/helper01_SCB-PCA-maker.R")


### Differential Expression Analysis on Sitka Spruce Cell type 
### experiment with limma + voom

# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


# Import file with tximport
quant.files <-
  dir("data/SCB-Salmon-Oct19",
      pattern = "quant.sf",
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE)

samples <-
  quant.files %>%
  str_replace("data/SCB-Salmon-Oct19/", "") %>%
  str_replace("_novRun_salmon_quant\\.sf", "")

# Rename samples.  
# 898 = R (resistance genotype)
# 903 = S (susceptible genotype)
samples <- samples %>% str_replace("898", "R")
samples <- samples %>% str_replace("903", "S")

samples_2_files <-
  data.frame(samples, quant.files,
             stringsAsFactors = FALSE)

txi <-
  tximport(samples_2_files[, "quant.files"],
           type = "salmon", txOut = TRUE, varReduce = TRUE,
           countsFromAbundance = "lengthScaledTPM")

colnames(txi$counts) <- samples_2_files$samples



# Create experimental design
expDes <- 
  data.frame(
    sample = samples,
    gType = factor(c(samples %>% str_extract("^\\w")),
                   levels = c("S", "R")),
    cType = factor(c(samples %>% 
                       str_extract("_\\w{2,3}_") %>% 
                       str_replace_all("_", "")),
                            levels = c("CP", "DSC")),
             bioRep = as.numeric(samples %>% str_extract("\\d$")))

expDes$group <- with(expDes, interaction(gType, cType))

#        sample gType cType bioRep group
# 1   R_CP_rep1     R    CP      1  R.CP
# 2   R_CP_rep2     R    CP      2  R.CP
# 3   R_CP_rep3     R    CP      3  R.CP
# 4  R_DSC_rep1     R   DSC      1 R.DSC
# 5  R_DSC_rep2     R   DSC      2 R.DSC
# 6  R_DSC_rep3     R   DSC      3 R.DSC
# 7   S_CP_rep1     S    CP      1  S.CP
# 8   S_CP_rep2     S    CP      2  S.CP
# 9   S_CP_rep3     S    CP      3  S.CP
# 10 S_DSC_rep1     S   DSC      1 S.DSC
# 11 S_DSC_rep2     S   DSC      2 S.DSC
# 12 S_DSC_rep3     S   DSC      3 S.DSC


# Load counts into DGEList object from edgeR package.
# x <- DGEList(counts = raw, group = expDes$group)
x <- DGEList(counts = txi$counts, group = expDes$group)

# Keep only genes with at least 1 count-per-million reads (cpm) in 
# at least 3 samples
dim(x)
# [1] 101973     12

# Filter low-expression contigs
y <- x[(rowSums(cpm(x) > 1) >= 3), ]

dim(y)
# [1] 26787    12

# Reset depth
y$samples$lib.size <- colSums(y$counts)

# TMM Normalization by Depth
y <- calcNormFactors(y)

scb_cpm <- cpm(y)

write.table(
  scb_cpm,
  "results/SCB.normalized_cpm.lowExpFiltered.txt",
  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


# make model matrix
# Interaction design
modMat <- model.matrix(~ gType * cType, expDes)
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat) <- gsub(":", "_", colnames(modMat))

#    Intercept gTypeR cTypeDSC gTypeR_cTypeDSC
# 1          1      1        0               0
# 2          1      1        0               0
# 3          1      1        0               0
# 4          1      1        1               1
# 5          1      1        1               1
# 6          1      1        1               1
# 7          1      0        0               0
# 8          1      0        0               0
# 9          1      0        0               0
# 10         1      0        1               0
# 11         1      0        1               0
# 12         1      0        1               0


cont_matrix <- makeContrasts(
  # 1. DE between CP and DSC in susceptible genotype (Q903)
  S_cType = cTypeDSC,
  # 2. DE between CP and DSC in resistance genotype (H898)
  R_cType = cTypeDSC + gTypeR_cTypeDSC,
  # 3. DE in CP between genotype
  CP_gType = gTypeR,
  # 4. DE in DSC between genotype
  DSC_gType = gTypeR + gTypeR_cTypeDSC,
  levels = modMat)


# voom transformation
# voom-plot
# v <- voom(x, modMat, plot = FALSE)
v <- voom(y, modMat, plot = TRUE)


# p <- PCA_maker(expDes, v)
# ggsave("results/figures/Fig1a-SCB-PCA.30Jun.svg", plot = p, 
       # height = 6, width = 6)


# Linear modelling
fit <- lmFit(v, modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

summary(decideTests(fit3, method = "separate", adjust.method = "fdr", 
                    p.value = pCutoff, lfc = lfcCutoff))
#        S_cType R_cType CP_gType DSC_gType
# Down       424     404     1530      1538
# NotSig   26119   25520    23903     23286
# Up         244     863     1354      1963

focus_terms <- colnames(cont_matrix)

results <-
  map_df(focus_terms, function(f) {
    tmp <- topTable(fit3, coef = f, number = Inf, 
                    sort.by = "none", adjust.method = "fdr")
    tmp$focus_term <- f
    tmp$cds <- row.names(tmp)
    return(tmp)
  })

str(results)
# 'data.frame':	107088 obs. of  8 variables:

table(results$focus_term)
#  CP_gType DSC_gType   R_cType   S_cType 
#     26787     26787     26787     26787 

# Reorganize columns
results <- results %>% 
  select(cds, logFC, adj.P.Val, focus_term)

write.table(
  results, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
  "results/SCB.all_stats.18March.txt"
)


### Write out all DE contig ids with abs(logFC) >= 2 & adj. p-value <= 0.05
sigDE <- results %>% 
  filter(adj.P.Val <= pCutoff & abs(logFC) >= lfcCutoff)

write.table(
  sigDE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE,
  "results/SCB.sigDE_stats.18March.txt"
  )
