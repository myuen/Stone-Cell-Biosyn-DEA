library(plyr)

source("analysis/helper03-figure_maker.R")

results.long <- 
  read.delim("results/StoneCellBiosyn_pooledRun_stats.18aug.long.txt", header = TRUE)


lfcCutoff <- 2
pCutoff <- 0.01


summarizeRes <- function(results.long) {
  ddply(results.long, ~ focus, 
        function(x) {
          x <- subset(x, global.adj.P <= pCutoff);
          gt <- table(x$logFC >= lfcCutoff)["TRUE"]
          lt <- table(x$logFC < (-1 * lfcCutoff))["TRUE"]
          return(c("up" = gt, "down" = lt))
        })
}

summarizeRes(results.long)
#        focus up.TRUE down.TRUE
# 1   CP_gType    3705      3790
# 2  DSC_gType    4155      3733
# 3 H898_cType    1137       436
# 4 Q903_cType     438       519


### Create Venn diagram for contigs up-regulated in H898 and Q903
H898_DSC_upReg <- subset(results.long, results.long$focus == "H898_cType" &
                           results.long$logFC >= lfcCutoff & 
                           results.long$global.adj.P <= pCutoff)
dim(H898_DSC_upReg)
# [1] 1137    5

Q903_DSC_upReg <- subset(results.long, results.long$focus == "Q903_cType" &
                           results.long$logFC >= lfcCutoff & 
                           results.long$global.adj.P <= pCutoff)
dim(Q903_DSC_upReg)
# [1] 438   5







CP_gType_vp <- volcanoPlot(
  subset(results.long, results.long$focus == "CP_gType"), lfc, pCutoff, 
  "Volcano Plot for Pooled Assembly\nCortical Parenchyma Comparison Between Genotype")
ggsave(filename = "results/figures/CP_gtype_pooledRun_volcanoPlot.19aug.png",
  plot = CP_gType_vp, height = 8.5, width = 11, unit = "in")

DSC_gType_vp <- volcanoPlot(
  subset(results.long, results.long$focus == "DSC_gType"), lfc, pCutoff, 
  "Volcano Plot for Pooled Assembly\nDeveloping Stone Cells Comparison Between Genotype")
ggsave(filename = "results/figures/DSC_gtype_pooledRun_volcanoPlot.19aug.png",
       plot = DSC_gType_vp, height = 8.5, width = 11, unit = "in")

H898_cType_vp <- volcanoPlot(
  subset(results.long, results.long$focus == "H898_cType"), lfc, pCutoff, 
  "Volcano Plot for Pooled Assembly\nCelltype Comparison in H898")
ggsave(filename = "results/figures/H898_ctype_pooledRun_volcanoPlot.19aug.png",
       plot = H898_cType_vp, height = 8.5, width = 11, unit = "in")

Q903_cType_vp <- volcanoPlot(
  subset(results.long, results.long$focus == "Q903_cType"), lfc, pCutoff, 
  "Volcano Plot for Pooled Assembly\nCdelltype Comparison in Q903")
ggsave(filename = "results/figures/Q903_ctype_pooledRun_volcanoPlot.19aug.png",
       plot = Q903_cType_vp, height = 8.5, width = 11, unit = "in")


CP_gType_lfc <- plotFC2AvgExpr(
  subset(results.long, results.long$focus == "CP_gType"), lfc, pCutoff, 
  "LogFC Plot for Pooled Assembly\nCortical Parenchyma Comparison Between Genotype")
ggsave(filename = "results/figures/CP_gType_pooledRun_logFC.19aug.png",
       plot = CP_gType_lfc, height = 8.5, width = 11, unit = "in")

DSC_gType_lfc <- plotFC2AvgExpr(
  subset(results.long, results.long$focus == "DSC_gType"), lfc, pCutoff, 
  "LogFC Plot for Pooled Assembly\nDeveloping Stone Cells Comparison Between Genotype")
ggsave(filename = "results/figures/DSC_gType_pooledRun_logFC.19aug.png",
       plot = DSC_gType_lfc, height = 8.5, width = 11, unit = "in")

H898_cType_lfc <- plotFC2AvgExpr(
  subset(results.long, results.long$focus == "H898_cType"), lfc, pCutoff, 
  "LogFC Plot for Pooled Assembly\nCelltype Comparison in H898")
ggsave(filename = "results/figures/H898_cType_pooledRun_logFC.19aug.png",
       plot = H898_cType_lfc, height = 8.5, width = 11, unit = "in")

Q903_cType_lfc <- plotFC2AvgExpr(
  subset(results.long, results.long$focus == "Q903_cType"), lfc, pCutoff, 
  "LogFC Plot for Pooled Assembly\nCelltype Comparison in Q903")
ggsave(filename = "results/figures/Q903_cType_pooledRun_logFC.19aug.png",
       plot = Q903_cType_lfc, height = 8.5, width = 11, unit = "in")

