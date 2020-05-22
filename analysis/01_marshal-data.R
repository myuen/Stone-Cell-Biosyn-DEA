library(dplyr)
library(purrr)
library(stringr)


# Marshal output files from Salmon
rawFiles <- list.files("data/Salmon-Oct19", full.names = TRUE)

## Extract library namne from filename
libNames <- str_extract(basename(rawFiles), "\\d{3}_\\w{2,3}_\\w{4}")

## Replace genotype name.  898 = Resistance; 903 = Susceptible
libNames <- gsub("898", "R", libNames)
libNames <- gsub("903", "S", libNames)

## Use as names for good side effects later
names(rawFiles) <- libNames


## Read in all Salmon data
system.time(
  myBigList <- map(rawFiles, function(x) {
    content <- 
      read.table(x, header = TRUE, colClasses = c("character", rep("numeric", 4)))
    
    # Only extract columns that will be use for analysis
    content <- select(content, Name, NumReads)
    
    # Rename colname to "CDS" and respective library name
    colnames(content) <- c("CDS", names(x))
    return(content)
  })
)
#    user  system elapsed 
#   5.117   0.092   5.233 

str(myBigList)
# List of 12


# Flatten big list of 5 into single data frame
rawSalmonCounts <- flatten_dfr(myBigList)
colnames(rawSalmonCounts) <- c("CDS", names(myBigList))
str(rawSalmonCounts)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	101973 obs. of  13 variables:


write.table(rawSalmonCounts, "data/consolidated-Salmon-counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## write design matrix
expDes <- data.frame(sample = libNames,
                     gType = factor(c(rep("R", 6), rep("S", 6)),
                                    levels = c("S", "R")),
                     cType = factor(rep(c(rep("CP", 3), rep("DSC", 3)), 2),
                                    levels = c("CP", "DSC")),
                     bioRep = as.numeric(rep(c(1, 2, 3), 4)))
expDes$group <- with(expDes, interaction(gType, cType))
str(expDes)
expDes

write.table(expDes, "data/stone_cell_biosyn_exp_design.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
