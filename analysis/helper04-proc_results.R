# Takes a MArrayLM fitted model object, 
# 1 - generate "global" adjusted p-value, then
# 2 - run topTable on all coefficients, 
# 3 - glue adj. p-value to each toptable result set
# 4 - merge and return a single data frame

proc_results <- function(fit) {
  require(reshape2)
  
  # Generate "global" adjusted p-value
  global.adj.P <- fit$p.value
  global.adj.P[] <- p.adjust(global.adj.P, method = "fdr")
  global.adj.P <- as.data.frame(global.adj.P)
  global.adj.P$contig <- rownames(global.adj.P)
  global.adj.P <- melt(global.adj.P, id.vars = "contig", 
                       variable.name = "focus", value.name = "global.adj.P")
  
  allTopTableRes <- topTable(fit, number = Inf, adjust.method = "fdr", sort.by = "none")
  allTopTableRes$contig <- rownames(allTopTableRes)
  allTopTableRes <- melt(allTopTableRes, id.vars = c("contig", "AveExpr"), 
                         measure.vars = colnames(fit$contrasts), 
                         variable.name = "focus", value.name = "logFC")
  
  mergeRes <- merge(allTopTableRes, global.adj.P, 
                    by.x = c("contig", "focus"), by.y = c("contig", "focus"))
  mergeRes <- mergeRes[, c("contig", "logFC", "AveExpr", "global.adj.P", "focus")]  
    
}