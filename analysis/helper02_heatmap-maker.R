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
  
  # htmlwidgets::saveWidget(
  #   hm,
  #   file = paste("results/figures/", outfile, ".html", sep = ""),
  #   selfcontained = TRUE)
  
  return (hm)
}


# Takes a cpm in dataframe and output file as argument
# TimeCourse heatmap is base on dendrogram created from LMD data
make2HM <- function(cpm.1, cpm.2, out_prefix) {
  
  hm.dist <- as.dist(1- cor(t(cpm.1), method = "pearson"))
  hm.clust <- hclust(hm.dist, method = "complete")
  hm.dendrogram <- as.dendrogram(hm.clust)
  
  hm.1 <-
    heatmaply(cpm.1,
              colors = heat.colors(50),
              scale = "row",
              
              # Rowv = scb_dendrogram,
              Rowv = hm.dendrogram,
              show_dendrogram = c(TRUE, FALSE),
              # show_dendrogram = c(FALSE, FALSE),
              row_dend_left = TRUE,
              # row_dend_left = FALSE,
              
              hide_colorbar = TRUE,
              dendrogram = "row",
              showticklabels = c(TRUE, FALSE)
    )
  
  hm.2 <-
    heatmaply(cpm.2,
              colors = heat.colors(50),
              scale = "row",
              
              # Rowv = scb_dendrogram,
              Rowv = hm.dendrogram,
              show_dendrogram = c(FALSE, FALSE),
              # show_dendrogram = c(TRUE, FALSE),
              row_dend_left = FALSE,
              # row_dend_left = TRUE,
              
              hide_colorbar = TRUE,
              dendrogram = "row",
              
              showticklabels = c(TRUE, FALSE)
    )
  
  merged <- subplot(hm.1, hm.2, margin = 0)
  
  # htmlwidgets::saveWidget(merged,
  #                         file = paste("results/figures/", out_prefix, ".html", sep = ""),
  #                         # file = "results/figures/test.html",
  #                         selfcontained = TRUE)
  
  return (merged)
}
