library(dplyr)
library(purrr)
library(stringr)


# Marshal output files from Salmon
rawFiles <- list.files("data/Salmon-results/", full.names = TRUE)

## Extract library namne from filename
libNames <- str_extract(basename(rawFiles), "\\d{3}_\\w{2,3}_\\w{4}")

## Add prefix to library name
libNames <- gsub("898", "H898", libNames)
libNames <- gsub("903", "Q903", libNames)


## Use as names for good side effects later
names(rawFiles) <- libNames


## Read in all Salmon data
system.time(
  myBigList <- map(rawFiles, function(x) {
    content <- read.table(x, header = TRUE, colClasses = c("character", rep("numeric", 4)))
    
    # Only extract columns that will be use for analysis
    content <- select(content, Name, TPM)
    
    # Rename colname to "CDS" and respective library name
    colnames(content) <- c("CDS", names(x))
    return(content)
  })
)
#    user  system elapsed 
#   2.534   0.036   2.579 

  
str(myBigList)
# 
  
# Flatten big list of 5 into single data frame
rawSalmonCounts <- flatten_dfr(myBigList)
colnames(rawSalmonCounts) <- c("CDS", names(myBigList))
str(rawSalmonCounts)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':     87266 obs. of  6 variables:


write.table(rawSalmonCounts, "data/consolidated-Salmon-counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## write design matrix
expDes <- data.frame(sample = libNames,
                     gType = factor(c(rep("H898", 6), rep("Q903", 6)),
                                    levels = c("Q903", "H898")),
                     cType = factor(rep(c(rep("CP", 3), rep("DSC", 3)), 2),
                                    levels = c("CP", "DSC")),
                     bioRep = as.numeric(rep(c(1, 2, 3), 4)))
expDes$group <- with(expDes, interaction(gType, cType))
str(expDes)
expDes

write.table(expDes, "data/stone_cell_biosyn_exp_design.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
