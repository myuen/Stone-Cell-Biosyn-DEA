library(plyr)
source("analysis/helper05-figure_maker.R")

julyRes.long <- read.delim("results/StoneCellBiosyn_julyRun_DE_results.long.txt", header = TRUE)
julyRes.wide <- read.delim("results/StoneCellBiosyn_julyRun_DE_results.wide.txt", header = TRUE)

novRes.long <- read.delim("results/StoneCellBiosyn_novRun_DE_results.long.txt", header = TRUE)
novRes.wide <- read.delim("results/StoneCellBiosyn_novRun_DE_results.wide.txt", header = TRUE)

lfc <- 2
pCutoff <- 0.01

summarizeRes <- function(results.long) {
  ddply(results.long, ~ focus, 
        function(x) {
          x <- subset(x, global.adj.P < pCutoff);
          gt <- table(x$logFC > lfc)["TRUE"]
          lt <- table(x$logFC < (-1 * lfc))["TRUE"]
          return(c("up" = gt, "down" = lt))
        })
}

summarizeRes(julyRes.long)
#        focus up.TRUE down.TRUE
# 1   CP_gType     772       837
# 2  DSC_gType    1039       871
# 3 H898_cType     321       119
# 4 Q903_cType      76        77

summarizeRes(novRes.long)
#        focus up.TRUE down.TRUE
# 1   CP_gType     931      1127
# 2  DSC_gType    1062      1068
# 3 H898_cType     362       142
# 4 Q903_cType      88       163


jvp <- volcanoPlot(julyRes.long, lfc, pCutoff, "Volcano Plot for July Assembly")
ggsave(filename = paste("results/figures/julyRun_volcanoPlot.png", sep = ""), 
  plot = jvp, height = 8.5, width = 11, unit = "in")

jlfc <- plotFC2AvgExpr(julyRes.long, lfc, pCutoff, "July Assembly")
ggsave(filename = paste("results/figures/julyRun_logFC.png", sep = ""), 
       plot = jlfc, height = 8.5, width = 11, unit = "in")

nvp <- volcanoPlot(julyRes.long, lfc, pCutoff, "Volcano Plot for Nov Assembly")
ggsave(filename = paste("results/figures/novRun_volcanoPlot.png", sep = ""), 
       plot = jvp, height = 8.5, width = 11, unit = "in")

nlfc <- plotFC2AvgExpr(novRes.long, lfc, pCutoff, "Nov Assembly")
ggsave(filename = paste("results/figures/novRun_logFC.png", sep = ""), 
       plot = nlfc, height = 8.5, width = 11, unit = "in")
