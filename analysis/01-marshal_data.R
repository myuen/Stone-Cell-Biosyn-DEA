source("analysis/helper01-consolidate_salmon_results.R")


# Marshal output files from Sailfish
julyRunFileList <- list.files("data/julyRun", 
                              pattern = "*_julyRun.quant.sf", 
                              full.names = TRUE)

test_that("Exactly 12 Salmon/Sailfish output files from July assembly are found",
          expect_equal(12, length(julyRunFileList)))

july_rRNA_id <- scan("data/julyRun/putative_rRNA.id", what = "")

julyCounts <- consolidate_results(julyRunFileList, july_rRNA_id)

###

novRunFileList <- list.files("data/novRun", 
                              pattern = "*_novRun.quant.sf", 
                              full.names = TRUE)

test_that("Exactly 12 Salmon/Sailfish output files from assembly are found",
          expect_equal(12, length(novRunFileList)))

nov_rRNA_id <- scan("data/novRun/putative_rRNA.id", what = "")

novCounts <- consolidate_results(novRunFileList, nov_rRNA_id)


###

write.table(julyCounts, "data/julyRun/consolidated-julyRun-Salmon-results.txt",
            sep = "\t", quote = FALSE)

write.table(novCounts, "data/novRun/consolidated-novRun-Salmon-results.txt",
            sep = "\t", quote = FALSE)
