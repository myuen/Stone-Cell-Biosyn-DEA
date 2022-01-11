require(heatmaply)

makeHM <- function(hm.cpm, out_prefix) {
  
  hm.cpm <- as.matrix(hm.cpm)
  hm.dist <- as.dist(1- cor(t(hm.cpm), method = "pearson"))
  hm.clust <- hclust(hm.dist, method = "complete")
  hm.dendrogram <- as.dendrogram(hm.clust)
  
  hm <- heatmaply(hm.cpm,
                  colors = heat.colors(50),
                  scale = "row",
                  
                  Rowv = hm.dendrogram,
                  show_dendrogram = c(TRUE, FALSE),
                  row_dend_left = TRUE,
                  
                  hide_colorbar = TRUE,
                  dendrogram = "row",
                  showticklabels = c(TRUE, FALSE)
                  # file = outfile
  )
  
  # orca(hm, file = outfile, scale = "none")
  
  htmlwidgets::saveWidget(
    hm,
    file = paste("results/figures/", out_prefix, ".html", sep = ""),
    selfcontained = TRUE)
  
  return (hm)
}
