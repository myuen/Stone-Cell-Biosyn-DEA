library(plyr)

source("analysis/helper03-figure_maker.R")

results.long <- 
  read.delim("results/StoneCellBiosyn_pooledRun_stats.18aug.long.txt", header = TRUE)


lfc <- 2
pCutoff <- 0.01


summarizeRes <- function(results.long) {
  ddply(results.long, ~ focus, 
        function(x) {
          x <- subset(x, global.adj.P <= pCutoff);
          gt <- table(x$logFC >= lfc)["TRUE"]
          lt <- table(x$logFC < (-1 * lfc))["TRUE"]
          return(c("up" = gt, "down" = lt))
        })
}

summarizeRes(results.long)
#        focus up.TRUE down.TRUE
# 1   CP_gType     772       837
# 2  DSC_gType    1039       871
# 3 H898_cType     321       119
# 4 Q903_cType      76        77


pvp <- volcanoPlot(results.long, lfc, pCutoff, "Volcano Plot for Pooled Assembly")
ggsave(filename = "results/figures/pooledRun_volcanoPlot.18aug.png",
  plot = pvp, height = 8.5, width = 11, unit = "in")

plfc <- plotFC2AvgExpr(results.long, lfc, pCutoff, "Pooled Assembly")
ggsave(filename = "results/figures/pooledRun_logFC.18aug.png",
       plot = plfc, height = 8.5, width = 11, unit = "in")
