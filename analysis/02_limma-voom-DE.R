require(edgeR)
require(ggplot2)
require(limma)
require(tibble)


pCutoff <- 0.01
lfcCutoff <- 2


# source("analysis/helper02-DE_workflow.R")
source("analysis/helper01-PCA_maker.R")
source("analysis/helper02-proc_results.R")


# Load raw Salmon quantification results
raw <- read.table("data/consolidated-Salmon-counts.txt", header = TRUE, 
                  colClasses = c("character", rep("numeric", 12)))
raw <- column_to_rownames(raw, "CDS")
str(raw)
# 'data.frame':	101973 obs. of  12 variables:


# Load experimental design
expDes <- read.table("data/stone_cell_biosyn_exp_design.tsv", header = TRUE)
expDes$gType <- relevel(expDes$gType, ref = "Q903")


# Load counts into DGEList object from edgeR package.
x <- DGEList(counts = raw, group = expDes$group)


# Keep only genes with at least 1 count-per-million reads (cpm) in at least 3 samples
dim(x)
# [1] 101973     12


# Filter low-expression contigs
x <- x[(rowSums(cpm(x) > 1) >= 3), ]


dim(x)
# [1] 26806    12


# Reset depth
x$samples$lib.size <- colSums(x$counts)


# TMM Normalization by Depth
x <- calcNormFactors(x)


write.table(cpm(x), "results/StoneCellBiosyn_pooledRun.normalized_cpm.lowExpFiltered.txt",
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


# MDS analysis
png("results/figures/pooledRun_MDS.25oct.png", height = 1200, width = 1680)
plotMDS(x, top = Inf)
dev.off()


# make model matrix
# Interaction design
modMat <- model.matrix(~ gType * cType, expDes)
colnames(modMat)
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat)
cont_matrix <-
  makeContrasts(
    CP_gType = gTypeH898,
    DSC_gType = gTypeH898 + gTypeH898_cTypeDSC, # H898 DSC - Q903 DSC
    Q903_cType = cTypeDSC,
    H898_cType = cTypeDSC + gTypeH898_cTypeDSC,
    levels = modMat)


# voom transformation
# voom-plot
v <- voom(x, modMat, plot = TRUE)


p <- PCA_maker(expDes, v)
ggsave("results/figures/pooledRun_PCA.25oct.png",
       plot = p, height = 8.5, width = 11, unit = "in")


# Linear modelling
fit <- lmFit(v, modMat)

fit <- contrasts.fit(fit, cont_matrix)

fit <- eBayes(fit)


# N.B.  LogFC cut off threashold set at 2
summary(decideTests(fit, method = "global", adjust.method = "fdr", 
                    p.value = pCutoff, lfc = lfcCutoff))
#        CP_gType DSC_gType Q903_cType H898_cType
# Down       1223      1200        348        276
# NotSig    24614     24078      26245      25722
# Up          969      1528        213        808

results <- proc_results(fit)


summary(results$focus)
#   CP_gType  DSC_gType Q903_cType H898_cType 
#      39839      39839      39839      39839 


results.wide <- reshape(results, direction = "wide", 
                        idvar = "contig", timevar = "focus")

### Write results files
write.table(
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
  results, "results/StoneCellBiosyn_pooledRun_stats.25oct.long.txt"
)

write.table(
  results.wide, quote = FALSE,  sep = "\t", row.names = FALSE,
  "results/StoneCellBiosyn_pooledRun_stats.25oct.wide.txt"
)


### Write out all DE contig ids with abs(logFC) >= 2 & adj. p-value <= 0.01
write.table(
  subset(results, results$global.adj.P <= pCutoff & abs(results$logFC) >= lfcCutoff),
  "results/StoneCellBiosyn_pooledRun.sigDE.25oct.txt", quote = FALSE,
  sep = "\t", col.names = TRUE, row.names = FALSE)
