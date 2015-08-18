## function to make Volcano plot
volcanoPlot <- function(df, logFC, pCutoff, title) {
  require (ggplot2)
  
  p <- ggplot(data = subset(df, (abs(logFC) < lfc | global.adj.P > pCutoff)), 
         aes(x = logFC, y = -log10(global.adj.P))) + geom_point(colour = "black", size = 1) + 
    geom_point(data = subset(df, (abs(logFC) >= lfc & global.adj.P <= pCutoff)), 
               aes(logFC, -log10(global.adj.P)), color = "red1", size = 1) + 
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "red") + 
    geom_vline(xintercept = lfc, colour = "red") + 
    geom_vline(xintercept = -(lfc), colour = "red") +
    labs(title = title, x = "Log 2 Fold Change", y = "-log 10 (Adjusted P Value)") + 
    scale_x_continuous(limits = c(-20, 20)) + theme_bw() + 
    theme(plot.title = element_text(size = rel(2)), legend.position = "none")
  
  return (p)
}


## function to plot logFC against Average Expression
plotFC2AvgExpr <- function(df, lfc, pCutoff, title) {
  require (ggplot2)
  
  ggplot(subset(df, global.adj.P > pCutoff), aes(x = AveExpr, y = logFC)) + 
    # stat. insignificant data points
    geom_point(shape = 20, color = "#666666", size = 1) + 
    # stat. sig and up regulated data points
    geom_point(data = subset(df, global.adj.P <= pCutoff & logFC >= lfc),
               aes(x = AveExpr, y = logFC), 
               colour = "red", fill= "red", shape = 24, size = 0.8) + 
    # stat. sig and down regulated data points
    geom_point(data = subset(df, global.adj.P <= pCutoff & logFC <= (-1 * lfc)),
               aes(x = AveExpr, y = logFC), 
               colour = "blue", fill= "blue", shape = 25, size = 0.8) + 
    geom_abline(aes(intercept = 0, slope = 0), colour = "black") + 
    scale_y_continuous(limits = c(-20, 20)) + theme_bw() +
    labs(title = title, x = "Log 2 Average Expression Counts", y = "Log 2 Fold Changes") +
    theme(plot.title = element_text(size = rel(2)))
}
