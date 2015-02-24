consolidate_results <- function(fileList, rRNA) {
  require(testthat)
  require(plyr)

  test_that("Exactly 12 Salmon/Sailfish output files from assembly are found",
            expect_equal(12, length(fileList)))
  
  ## Extract library names from list of file names to be used in data file
  libNames <- gsub("data/[a-z]+Run/", "", fileList)
  libNames <- gsub("_[a-z]+Run.quant.sf", "", libNames)
  libNames <- gsub("898", "H898", libNames)
  libNames <- gsub("903", "Q903", libNames)
  str(libNames)
  

  ## Read one file to learn how many rows we expect and to grab rownames
  tmp <- read.table(fileList[1], row.names = 1,
                    ## specifying colClasses speeds this up 2x
                    colClasses = rep(c("character", "numeric"), c(1, 4)))
  (n <- nrow(tmp))
  contigIds <- rownames(tmp)
  

  ## Read in all Sailfish data
  system.time(rawSalmonCounts <- aaply(fileList, 1, function(x) {
    mDat <- 
      read.table(x, row.names = 1, nrows = n,
                 ## specifying colClasses speeds this up 2x
                 colClasses = rep(c("character", "numeric"), c(1, 4)))
    return(mDat$V5)})
  )

  rawSalmonCounts <- data.frame(t(rawSalmonCounts))
  colnames(rawSalmonCounts) <- libNames
  rownames(rawSalmonCounts) <- contigIds
  str(rawSalmonCounts)
  
  summary(rownames(rawSalmonCounts) %in% rRNA)
  # Mode   FALSE    TRUE    NA's 
  # logical  369946       8       0
  rawSalmonCounts <- rawSalmonCounts[!(rownames(rawSalmonCounts) %in% rRNA),]
  summary(rownames(rawSalmonCounts) %in% rRNA)
  # Mode   FALSE    NA's 
  # logical  369946       0

  return(rawSalmonCounts) 
}
