PCA_maker <- function(expDes, v) {
  require(ggplot2)
    
  pca <- as.data.frame(prcomp(t(v$E))$x)
  pca$sample <- gsub("_salmon", "", rownames(pca))
  pca <- merge(pca, expDes, by.x = "sample", by.y = "sample")
  
  
  p <- ggplot(data = pca, aes_string(x = "PC1", y = "PC2", color = "gType",  shape = "cType")) + 
    geom_point(size = 4) +
    geom_text(data = pca, aes(x = PC1, y = PC2, label = sample), 
              size = 3, vjust = 2.25, hjust = 0.5) + 
    theme_bw() + theme(legend.background = element_rect(size = 5, linetype = "solid"),
                       legend.justification = c(0, 0),
                       legend.position = c(0.85, 0.05),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12)) +
    scale_color_manual(values = c("Q903" = "#e34a33", "H898" = "#2b8cbe")) +
    labs(title=paste("Principal Component Analysis (PCA) of Stone Cell Biosynthesis 
                   RNA-Seq Libraries from"),
         shape = "Celltype", color = "Genotype", size = 40)
  
  return(p)
}