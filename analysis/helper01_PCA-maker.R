PCA_maker <- function(expDes, v) {
  require(ggplot2)

  pca <- as.data.frame(prcomp(t(v$E))$x)
  pca <- merge(pca, expDes, by.x = 0, by.y = "sample")
  
  p <- ggplot(
    data = pca, 
    aes_string(x = "PC1", y = "PC2", color = "group",  shape = "group")) +
    geom_point(size = 4) +

    scale_color_manual(values = c("S.CP" = "#e34a33",
                                  "S.DSC" = "#e34a33",
                                  "R.CP" = "#2b8cbe",
                                  "R.DSC" = "#2b8cbe")) +
    scale_shape_manual(values = c("S.CP" = 16, 
                                  "S.DSC" = 17,
                                  "R.CP" = 16, 
                                  "R.DSC" = 17)) +

    # Text on labels
    geom_text(data = pca, aes(x = PC1, y = PC2, label = colnames(v$E)),
              size = 3, vjust = 2.25, hjust = 0.5) +
    theme_bw()

  return(p)
}
