library(dplyr)
library(purrr)
library(stringr)


# Marshal output files from Salmon for TimeCourse libraries
rawFiles <- list.files("data/TimeCourse-Salmon-7Apr2020", full.names = TRUE)

# Extract library namne from filename
libNames <- str_replace(basename(rawFiles), "_quant.sf", "")

# Add prefix to 898 (R) and 903 (S)
libNames <- str_replace(libNames, "898", "R")
libNames <- str_replace(libNames, "903", "S")

# Four timepoints: A6 = T1, A7 = T2, A8 = T3, A10 = T4
libNames <- str_replace(libNames, "A6", "T1")
libNames <- str_replace(libNames, "A7", "T2")
libNames <- str_replace(libNames, "A8", "T3")
libNames <- str_replace(libNames, "A10", "T4")


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
#   5.729   0.349   6.255 

rawSalmonCounts <- myBigList %>% 
  select(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24,
         26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48)

colnames(rawSalmonCounts) <- c("CDS", names(rawFiles))

# Re-organize column by timepoints
rawSalmonCounts <- rawSalmonCounts %>% select(
  "CDS", 
  "R_T1-1", "R_T1-2", "R_T1-3", "R_T2-1", "R_T2-2", "R_T2-3", 
  "R_T3-1", "R_T3-3", "R_T3-4", "R_T4-1", "R_T4-3", "R_T4-4", 
  "S_T1-1", "S_T1-2", "S_T1-4", "S_T2-4", "S_T2-7", "S_T2-8",
  "S_T3-1", "S_T3-2", "S_T3-3", "S_T4-1", "S_T4-3", "S_T4-4")

str(rawSalmonCounts)
# 'data.frame':	101973 obs. of  25 variables:

write.table(rawSalmonCounts, "data/SCTC-consolidated-Salmon-counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
