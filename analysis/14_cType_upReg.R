library(dplyr)

sigDE <- read.table("results/SCB.sigDE_annotated.15Jun.txt",
                    header = TRUE, sep = "\t", quote = "", 
                    stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	8318 obs. of  7 variables:


sigDE.annot <- sigDE %>% select(cds, accs, annots, pfam_accs)
sigDE.annot <- unique(sigDE.annot)
str(sigDE.annot)
# 'data.frame':	5307 obs. of  4 variables:


# Identify contigs that are up-regulated in S_cType and R_cType
R_cType_upReg <- sigDE %>% 
  filter(logFC > 0 & sigDE$focus == "R_cType")
str(R_cType_upReg)
# 'data.frame':	870 obs. of  7 variables:

S_cType_upReg <- sigDE %>% 
  filter(logFC > 0 & sigDE$focus == "S_cType")
str(S_cType_upReg)
# 'data.frame':	247 obs. of  7 variables:

# Identify contigs that are up-regulated in both genotype
interx <- intersect(R_cType_upReg$cds, S_cType_upReg$cds)
length(interx)
# [1] 75

R_excl <- setdiff(R_cType_upReg$cds, interx)
length(R_excl)
# [1] 795

S_excl <- setdiff(S_cType_upReg$cds, interx)
length(S_excl)
# [1] 172

cType_upReg <- 
  data.frame("cds" = R_excl, "set" = "R_excl", stringsAsFactors = FALSE)
cType_upReg <- 
  rbind(cType_upReg, data.frame("cds" = interx, "set" = "Intersection"))
cType_upReg <- 
  rbind(cType_upReg, data.frame("cds" = S_excl, "set" = "S_excl"))

str(cType_upReg)
# 'data.frame':	1042 obs. of  2 variables:

table(cType_upReg$set)
# Intersection       R_excl       S_excl 
#           75          795          172 


### Annotate cType up-regulated contigs with BLAST results
cType_upReg_annot <- left_join(cType_upReg, sigDE.annot)
# Joining, by = "cds"

str(cType_upReg_annot)
# 'data.frame':	1042 obs. of  5 variables:

write.table(cType_upReg_annot, 
            "results/SCB.upReg_DSC.annot.15Jun.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
