library(dplyr)
library(stringr)
library(tidyr)

library(org.At.tair.db)
# columns(org.At.tair.db)

library(GO.db)
# columns(GO.db)


### Read BLAST results for all DE contigs
blast <- read.table("results/SCB.upRegDSC.17Jun.blastpTAIR10.txt", 
                    header = TRUE, sep = "\t", 
                    quote = "", stringsAsFactors = FALSE)

blast <- blast %>% 
  dplyr::select(qseqid, sseqid) %>%
  distinct()

# Rename column for easier joining downstream
colnames(blast) <- c("cds", "TAIR")

# Strip TAIR10 version ID from accession number
blast$TAIR <- 
  blast$TAIR %>% 
  str_replace(".\\d$", "")

str(blast)
# 'data.frame':	1025 obs. of  3 variables:


### Map TAIR ID to GOID.  Filter for "Biological Process"
tair2go <- 
  select(org.At.tair.db, 
         keys = unique(blast$TAIR), 
         columns = c("GO", "ONTOLOGY")) %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::select(-c(EVIDENCE, ONTOLOGY))
# 'select()' returned 1:many mapping between keys and columns

str(tair2go)
# 'data.frame':	2436 obs. of  2 variables:


### Map GOID to GO terms
goid2term <- select(GO.db, 
                    keys = unique(tair2go$GO), 
                    columns = c("GOID", "TERM"))
# 'select()' returned 1:1 mapping between keys and columns


### Join table.  Map TAIR ID to GO terms
tair2goterm <- left_join(tair2go, goid2term, by = c("GO"= "GOID"))

str(tair2goterm)
# 'data.frame':	2436 obs. of  3 variables:

blast2goterm <- inner_join(blast, tair2goterm)
# Joining, by = "TAIR"

str(blast2goterm)
# 'data.frame':	4989 obs. of  4 variables:


sctc.s_peakedT2 <- 
  scan("results/SCTC.S_PeakedT2.txt", what = "character")
sctc.s_peakedT3 <- 
  scan("results/SCTC.S_PeakedT3.txt", what = "character")
sctc.r_peakedT2 <- 
  scan("results/SCTC.R_PeakedT2.txt", what = "character")
sctc.r_peakedT3 <- 
  scan("results/SCTC.R_PeakedT3.txt", what = "character")


### Group annotations
# blast.grouped <- 
#   blast2goterm %>% 
#   group_by(TERM)


### S T2
# sctc.s_peakedT2.go_terms <-
#   blast.grouped %>%
#   filter(cds %in% sctc.s_peakedT2) %>%
#   count(sort = TRUE, name = "count")

# str(sctc.s_peakedT2.go_terms)
# grouped_df [53 × 2] (S3: grouped_df/tbl_df/tbl/data.frame)

sctc.s_peakedT2.go_terms <-
  blast2goterm %>%
  filter(cds %in% sctc.s_peakedT2)

str(sctc.s_peakedT2.go_terms)
# 'data.frame':	100 obs. of  4 variables:

write.table(sctc.s_peakedT2.go_terms, 
            "results/SCTC.S_PeakedT2.GoTerms.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

### S T3
# sctc.s_peakedT3.go_terms <-
#   blast.grouped %>%
#   filter(cds %in% sctc.s_peakedT3) %>%
#   count(sort = TRUE, name = "count")

# str(sctc.s_peakedT3.go_terms)
# grouped_df [71 × 2] (S3: grouped_df/tbl_df/tbl/data.frame)

sctc.s_peakedT3.go_terms <-
  blast2goterm %>%
  filter(cds %in% sctc.s_peakedT3)

str(sctc.s_peakedT3.go_terms)
# 'data.frame':	96 obs. of  4 variables:

write.table(sctc.s_peakedT3.go_terms, 
            "results/SCTC.S_PeakedT3.GoTerms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T2
# sctc.r_peakedT2.go_terms <-
#   blast.grouped %>%
#   filter(cds %in% sctc.r_peakedT2) %>%
#   count(sort = TRUE, name = "count")
# 
# str(sctc.r_peakedT2.go_terms)
# grouped_df [239 × 2] (S3: grouped_df/tbl_df/tbl/data.frame)

sctc.r_peakedT2.go_terms <-
  blast2goterm %>%
  filter(cds %in% sctc.r_peakedT2)

str(sctc.r_peakedT2.go_terms)
# 'data.frame':	643 obs. of  4 variables:

write.table(sctc.r_peakedT2.go_terms, 
            "results/SCTC.R_PeakedT2.GoTerms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### R T3
# sctc.r_peakedT3.go_terms <-
#   blast.grouped %>%
#   filter(cds %in% sctc.r_peakedT3) %>%
#   count(sort = TRUE, name = "count")
# 
# str(sctc.r_peakedT3.go_terms)
# grouped_df [295 × 2] (S3: grouped_df/tbl_df/tbl/data.frame)

sctc.r_peakedT3.go_terms <-
  blast2goterm %>%
  filter(cds %in% sctc.r_peakedT3)

str(sctc.r_peakedT3.go_terms)
# 'data.frame':	1155 obs. of  4 variables:

write.table(sctc.r_peakedT3.go_terms, 
            "results/SCTC.R_PeakedT3.GoTerms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
