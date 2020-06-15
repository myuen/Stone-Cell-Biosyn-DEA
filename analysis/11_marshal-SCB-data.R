library(dplyr)
library(purrr)
library(stringr)


# Marshal output files from Salmon for LMD libraries
rawFiles <- list.files("data/LMD-Salmon-Oct19", full.names = TRUE)

## Extract library namne from filename
libNames <- str_extract(basename(rawFiles), "\\d{3}_\\w{2,3}_\\w{4}")

## Replace genotype name.  898 = Resistance; 903 = Susceptible
libNames <- gsub("898", "R", libNames)
libNames <- gsub("903", "S", libNames)

## Use as names for good side effects later
names(rawFiles) <- libNames


## Read in all Salmon data
system.time(
  myBigList <- map_dfc(rawFiles, function(x) {
    content <- 
      read.table(x, header = TRUE, colClasses = c("character", rep("numeric", 4)))
    
    # Only extract columns that will be use for analysis
    content <- content %>% select(Name, NumReads)
    
    # Rename colname to "CDS" and respective library name
    colnames(content) <- c("CDS", names(x))
    return(content)
  }))
#    user  system elapsed 
#   5.117   0.092   5.233 

rawSalmonCounts <- myBigList %>% 
  select(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)

colnames(rawSalmonCounts) <- c("CDS", names(rawFiles))

str(rawSalmonCounts)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	101973 obs. of  13 variables:

write.table(rawSalmonCounts, "data/SCB-consolidated-Salmon-counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
