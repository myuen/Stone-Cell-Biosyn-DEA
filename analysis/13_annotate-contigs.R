library(dplyr)
library(stringr)
library(tidyr)

library(org.At.tair.db)
# columns(org.At.tair.db)

library(GO.db)
# columns(GO.db)


### Read DEA results
sigDE <- read.table("results/SCB.sigDE_stats.17Jun.txt",
                    header = TRUE, stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	2283 obs. of  4 variables:


# Read file with contig ID up-regulated in developing stone cells (DSC)
upReg <- sigDE %>%
  filter(logFC >= 2) %>%
  dplyr::select(cds, focus_term)

str(upReg)
# 'data.frame':	1373 obs. of  2 variables:


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


# Join up-regulated contig ID with TAIR10 BLAST annotations
upReg_annot <- inner_join(upReg, blast)
# Joining, by = "cds"

str(upReg_annot)
# 'data.frame':	1086 obs. of  4 variables:


### Map TAIR ID to GOID.  Filter for "Biological Process"
tair2go <- 
  select(org.At.tair.db, 
         keys = unique(upReg_annot$TAIR), 
         columns = c("GO", "ONTOLOGY")) %>%
  dplyr::filter(ONTOLOGY == "BP") %>%
  dplyr::select(-c(EVIDENCE, ONTOLOGY))

str(tair2go)
# 'data.frame':	2436 obs. of  2 variables:


### Map GOID to GO terms
goid2term <- select(GO.db, 
                    keys = unique(tair2go$GO), 
                    columns = c("GOID", "TERM"))


### Join table.  Map TAIR ID to GO terms
tair2goterm <- left_join(tair2go, goid2term, by = c("GO"= "GOID"))

str(tair2goterm)
# 'data.frame':	2436 obs. of  3 variables:


### Join table.  Map up-regulated contigs to GO terms
upReg_annot <- inner_join(upReg_annot, tair2goterm)
# Joining, by = "TAIR"

str(upReg_annot)
# 'data.frame':	5157 obs. of  5 variables:

# write.table(upReg_annot, "results/upReg.go_terms.txt", 
#             quote = FALSE, row.names = FALSE, sep = "\t")


### Group annotations
upReg_annot.grouped <- 
  upReg_annot %>% 
  group_by(focus_term, TERM)
  

### Subset annotations for up-regulated contigs in susceptible genotype
s.go_terms <-
  upReg_annot.grouped %>%
  count(sort = TRUE, name = "count") %>%
  filter(focus_term == "S_cType")

str(s.go_terms)
# grouped_df [335 × 3] (S3: grouped_df/tbl_df/tbl/data.frame)

write.table(s.go_terms, "results/s.go_terms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### Subset annotations for up-regulated contigs in resistance genotype
r.go_terms <- 
  upReg_annot.grouped %>%
  count(sort = TRUE, name = "count") %>% 
  filter(focus_term == "R_cType")

str(r.go_terms)
# grouped_df [583 × 3] (S3: grouped_df/tbl_df/tbl/data.frame)

write.table(r.go_terms, "results/r.go_terms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


### Subset annotations for up-regulated contigs in both genotype
rs <- upReg[duplicated(upReg$cds), "cds"]

str(rs)
# chr [1:80]

rs.go_terms <-
  upReg_annot %>% 
  filter(cds %in% rs) %>% 
  group_by(TERM) %>% 
  count(sort = TRUE, name = "count")

write.table(rs.go_terms, "results/rs.go_terms.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
