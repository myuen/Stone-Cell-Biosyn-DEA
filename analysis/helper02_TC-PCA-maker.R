PCA_maker <- function(expDes, v) {
  require(ggplot2)
  
  pca <- as.data.frame(prcomp(t(v$E))$x)
  pca <- merge(pca, expDes, by.x = 0, by.y = "sample")
  
  p <- ggplot(
    data = pca, 
    aes_string(x = "PC1", y = "PC2", color = "group",  shape = "group")) +
    geom_point(size = 4) +
    
    scale_color_manual(values = c("S.T1" = "#e34a33",
                                  "S.T2" = "#e34a33",
                                  "S.T3" = "#e34a33",
                                  "S.T4" = "#e34a33",
                                  "R.T1" = "#2b8cbe",
                                  "R.T2" = "#2b8cbe",
                                  "R.T3" = "#2b8cbe",
                                  "R.T4" = "#2b8cbe")) +
    scale_shape_manual(values = c("S.T1" = 15, 
                                  "S.T2" = 16, 
                                  "S.T3" = 17, 
                                  "S.T4" = 18, 
                                  "R.T1" = 15,
                                  "R.T2" = 16,
                                  "R.T3" = 17, 
                                  "R.T4" = 18)) +
    
    # Text on labels
    geom_text(data = pca, aes(x = PC1, y = PC2, label = colnames(v$E)),
              size = 3, vjust = 2.25, hjust = 0.5) +
    theme_bw()
  
  return(p)
}
