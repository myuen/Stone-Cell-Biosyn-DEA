PCA_maker <- function(expDes, v) {
  require(ggplot2)

  pca <- as.data.frame(prcomp(t(v$E))$x)
  pca <- merge(pca, expDes, by.x = 0, by.y = "sample")
  
  p <- ggplot(
    data = pca, 
    aes_string(x = "PC1", y = "PC2", color = "gType",  shape = "cType"), label = colnames(v$E)) + 
    geom_point(size = 4) +
    # Text on labels
    geom_text(data = pca, aes(x = PC1, y = PC2, label = colnames(v$E)),
              size = 3, vjust = 2.25, hjust = 0.5) +
    theme_bw() + 
    scale_color_manual(values = c("S" = "#e34a33", "R" = "#2b8cbe"))

  return(p)
}
