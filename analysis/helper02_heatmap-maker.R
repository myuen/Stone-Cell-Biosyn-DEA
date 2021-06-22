require(heatmaply)

SCB_HM <- function(scb_cpm, sctc_cpm, outfile1, outfile2) {

  scb_cpm <- as.matrix(scb_cpm)
  
  scb_d <- as.dist(1- cor(t(scb_cpm), method = "pearson"))
  
  scb_c <- hclust(scb_d, method = "complete")
  
  scb_dendrogram <- as.dendrogram(scb_c)
  
  scb_hm <- 
    heatmaply(scb_cpm,
              colors = heat.colors(50),
              scale = "row",
              row_dend_left = TRUE,
              Rowv = scb_dendrogram,
              hide_colorbar = TRUE,
              dendrogram = "row",
              showticklabels = c(TRUE, FALSE),
              file = outfile1
    )
  
  sctc_cpm <- as.matrix(sctc_cpm)
  
  sctc_hm <- 
    heatmaply(sctc_cpm,
              colors = heat.colors(50),
              scale = "row",
              Rowv = scb_dendrogram,
              hide_colorbar = TRUE,
              dendrogram = "row",
              show_dendrogram = c(FALSE, FALSE),
              showticklabels = c(TRUE, FALSE),
              file = outfile2
    )
}