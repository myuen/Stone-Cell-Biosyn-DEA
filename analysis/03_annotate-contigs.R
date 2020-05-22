library(dplyr)
library(purrr)
library(rhmmer)
library(tidyr)


### Read DEA results
sigDE <- read.table("results/StoneCellBiosyn.sigDE_stats.21May.txt", 
                    header = TRUE, stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	8141 obs. of  4 variables:


### Read BLAST results for all DE contigs
blast <- read.table("results/sigDE.blastpNR.txt", header = TRUE, 
                    sep = "\t", quote = "", stringsAsFactors = FALSE)

blast <- blast %>% select(qseqid, sseqid, evalue, salltitles)

# Rename column for easier joining downstream
colnames(blast)[1] <- "cds"
str(blast)
# 'data.frame':	41437 obs. of  4 variables:


blast_flattened <- blast %>% 
  group_by(cds) %>% 
  nest()

blast_flattened$data <- map(blast_flattened$data, function(x){
  # Concatenate BLAST hits and delimit by ';' 
  x %>% select (sseqid, salltitles) %>% 
    summarise(accs = paste(sseqid, collapse=";"),
              annots = paste(salltitles, collapse = ";"))
})

blast_flattened <- blast_flattened %>% unnest(c(data))
str(blast_flattened)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	4347 obs. of  3 variables:


sigDE_annot <- left_join(sigDE, blast_flattened)
# Joining, by = "cds"

str(sigDE_annot)
# 'data.frame':	8141 obs. of  6 variables:


### Read HMMscan output
pfam <- read_domtblout("data/pfam.domtblout")
str(pfam)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':     857084 obs. of  23 variables:

# Only keep pfam with e-value less than 1e-10
pfam <- pfam %>% filter(sequence_evalue <= 1e-10) %>% 
  select(query_name, domain_accession) %>% unique()
str(pfam)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':     77458 obs. of  2 variables:

# Flatten pfam accession from multiple to a single line for each query
pfam_flattened <- pfam %>% 
  select (query_name, domain_accession) %>% 
  group_by(query_name) %>%
  summarise(pfam_accs = paste(domain_accession, collapse=";")) %>%
  as.data.frame()

# TransDecoder ID was renamed after HMMscan so we need 
# the conversion from old ID to new ID
id_transform <- read.delim("data/id_transform.txt", 
                           header = FALSE, stringsAsFactors = FALSE)
# Rename column for easier join later
colnames(id_transform) <- c("query_name", "cds")

pfam_flattened <- 
  left_join(pfam_flattened, id_transform) %>%
  select(cds, pfam_accs)
# Joining, by = "query_name"

str(pfam_flattened)
# 'data.frame': 43880 obs. of  2 variables:

sigDE_annot <- left_join(sigDE_annot, pfam_flattened)
# Joining, by = "cds"

str(sigDE_annot)
# 'data.frame':	8141 obs. of  7 variables:

# Supplementary Table 1
write.table(sigDE_annot, "results/StoneCellBiosyn.sigDE.annot.22May.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
