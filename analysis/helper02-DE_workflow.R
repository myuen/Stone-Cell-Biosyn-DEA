DE_workflow <- function(run, makePlots) {
  require(edgeR)
  require(limma)
  require(ggplot2)
  source("analysis/helper03-PCA_maker.R")
  
  raw <- read.table(paste("data/", run, "/consolidated-", run, "-Salmon-results.txt", 
                          sep = ""))
  
  n <- nrow(raw)
  
  if(run == "julyRun") {
    test_that("First Salmon output file has 369,954 rows",
              expect_equal(369946, n))
  } else {
    test_that("First Salmon output file has 104,132 rows",
              expect_equal(104132, n))
  }
  
  #' Load experimental design
  expDes <- read.table("data/stone_cell_biosyn_exp_design.tsv", header = TRUE)
  expDes$gType <- relevel(expDes$gType, ref = "Q903")
  
  
  #' Load counts into DGEList object from edgeR package.
  x <- DGEList(counts = raw, group = expDes$group)
  
  
  # Keep only genes with at least 1 count-per-million reads (cpm) in at least 4 samples
  dim(x)
  # for July assembly: [1] 369946     12
  # for Nov assembly: [1] 104132     12
  
  x <- x[(rowSums(cpm(x) > 1) >= 3), ]
  
  dim(x) 
  # for July assembly: [1] 35804    12
  # for Nov assembly: [1] 34407    12
  
  
  # Reset depth
  x$samples$lib.size <- colSums(x$counts)
  
  
  # TMM Normalization by Depth
  x <- calcNormFactors(x)
  
  
  # MDS analysis
  if(makePlots) {
    png(paste("results/figures/", run, "_MDS.png", sep = ""), height = 1200, width = 1680)
    plotMDS(x, top = Inf)
    dev.off()
  }
  
  
  # make model matrix
  # Interaction design
  modMat <- model.matrix(~ gType * cType, expDes)
  colnames(modMat)
  colnames(modMat) <- gsub("[()]", "", colnames(modMat))
  colnames(modMat) <- gsub(":", "_", colnames(modMat))
  colnames(modMat)
  cont_matrix <-
    makeContrasts(
      CP_gType = gTypeH898,
      DSC_gType = gTypeH898 + gTypeH898_cTypeDSC, # H898 DSC - Q903 DSC
      Q903_cType = cTypeDSC,
      H898_cType = cTypeDSC + gTypeH898_cTypeDSC,
      levels = modMat)
  
  
  #' voom transformation
  #+ voom-plot
  v <- voom(x, modMat, plot = TRUE)

  if (makePlots) {
    p <- PCA_maker(expDes, v)
    
    ggsave(
      filename = paste("results/figures/", run, "_PCA.png", sep = ""), 
      plot = p, height = 8.5, width = 11, unit = "in")
    
  }
  
  
  #' Linear modelling
  fit <- lmFit(v, modMat)
  fit <- contrasts.fit(fit, cont_matrix)
  fit <- eBayes(fit)
  
  return (fit)
}
