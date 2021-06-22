library(dplyr)
# library(purrr)
# library(rhmmer)
library(stringr)
library(tidyr)


### Read DEA results
sigDE <- read.table("results/SCB.sigDE_stats.17Jun.txt", 
                    header = TRUE, stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	2283 obs. of  4 variables:


# Identified contigs up-regulated in developing stone cells
upReg <- sigDE %>% 
  filter(logFC >= 2) %>%
  select(cds, focus_term)

str(upReg)
# 'data.frame':	1373 obs. of  2 variables:


### Read BLAST results for all DE contigs
blast <- read.table("results/SCB.upRegDSC.17Jun.blastpTAIR10.txt", 
                    header = TRUE, sep = "\t", 
                    quote = "", stringsAsFactors = FALSE)

blast <- blast %>% 
  select(qseqid, sseqid, salltitles) %>%
  distinct()

# Rename column for easier joining downstream
colnames(blast)[1:2] <- c("cds", "locus_name")

# Strip TAIR10 version ID from accession number
blast$locus_name <- 
  blast$locus_name %>% 
  str_replace(".\\d$", "")

str(blast)
# 'data.frame':	1025 obs. of  3 variables:


# Join up-regulated contig ID with TAIR10 BLAST annotations
sigDE_annot <- inner_join(upReg, blast)
# Joining, by = "cds"

str(sigDE_annot)
# 'data.frame':	1086 obs. of  4 variables:


# Refer to TAIR10 for column name header documentation
# https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO.README.txt
tair10_go <- read.delim("data/ATH_GO_GOSLIM.txt", 
                        sep = "\t", header = FALSE,  
                        comment.char = "!",
                        stringsAsFactors = FALSE)


# Column headers :explanation
# 1. locus name: standard AGI convention name
# 5. GO term: the actual string of letters corresponding to the GO ID
# 6. GO ID: the unique identifier for a GO term.  
# 8. Aspect: F=molecular function, C=cellular component, P=biological 13process. 
# 9. GOslim term: high level GO term helps in functional categorization.
tair10_go <- tair10_go[, c(1,5,6,8,9)]

colnames(tair10_go) <- 
  c("locus_name", "GO_term", "GO_id", "Aspect", "GO_slim")

tair10_go_bp <- tair10_go %>% filter(Aspect == "P")

str(tair10_go_bp)
# 'data.frame':	147211 obs. of  5 variables:

sigDE_annot <- inner_join(sigDE_annot, tair10_go_bp)
# Joining, by = "locus_name"

str(sigDE_annot)
# 'data.frame':	8944 obs. of  8 variables:

sigDE_annot %>% 
  group_by(focus_term, GO_term) %>% 
  select(GO_term) %>% 
  count(sort = TRUE) %>% 
  filter(focus_term == "S_cType")

sigDE_annot %>% 
  group_by(focus_term, GO_term) %>% 
  select(GO_term) %>% 
  count(sort = TRUE) %>% 
  filter(focus_term == "R_cType")

