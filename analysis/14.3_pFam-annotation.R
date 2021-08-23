library(dplyr)
library(rhmmer)

pfam <- read_domtblout("data/pfam.domtblout")

# Filter pfam by e-value
pfam <- pfam %>% 
  filter(sequence_evalue <= 1e-10)

str(pfam)
# tibble [204,061 × 23] (S3: tbl_df/tbl/data.frame)

# Remove pfam domain with no description
pfam <- pfam[!is.na(pfam$description),]

str(pfam)
# tibble [179,052 × 23] (S3: tbl_df/tbl/data.frame)


# TransDecoder ID was renamed.  Need to have a conversion 
# table from orig ID to renamed ID
origId_2_renamedId <- 
  read.csv("data/id_transform.txt", header = FALSE,
           sep = "\t", stringsAsFactors = FALSE)

colnames(origId_2_renamedId) <- c("orig", "renamed")


sctc.s_peakedT2 <- 
  scan("results/SCTC.S_PeakedT2.txt", what = "character")
# Read 31 items

sctc.s_peakedT2 <- 
  origId_2_renamedId[origId_2_renamedId$renamed %in% sctc.s_peakedT2,]

sctc.s_peakedT3 <- 
  scan("results/SCTC.S_PeakedT3.txt", what = "character")
# Read 26 items

sctc.s_peakedT3 <- 
  origId_2_renamedId[origId_2_renamedId$renamed %in% sctc.s_peakedT3,]


sctc.r_peakedT2 <- 
  scan("results/SCTC.R_PeakedT2.txt", what = "character")
# Read 158 items

sctc.r_peakedT2 <- 
  origId_2_renamedId[origId_2_renamedId$renamed %in% sctc.r_peakedT2,]


sctc.r_peakedT3 <- 
  scan("results/SCTC.R_PeakedT3.txt", what = "character")
# Read 342 items

sctc.r_peakedT3 <- 
  origId_2_renamedId[origId_2_renamedId$renamed %in% sctc.r_peakedT3,]


# sctc.s_peakedT2.pfam <- 
#   pfam %>% filter(query_name %in% sctc.s_peakedT2$orig) %>%
#   group_by(description) %>%
#   count(sort = TRUE, name = "count")

sctc.s_peakedT2.pfam <- pfam %>% 
  filter(query_name %in% sctc.s_peakedT2$orig) %>%
  select(query_name, domain_name, sequence_evalue, description)

str(sctc.s_peakedT2.pfam)
# tibble [24 × 4] (S3: tbl_df/tbl/data.frame)

write.table(sctc.s_peakedT2.pfam, 
            "results/SCTC.S_PeakedT2.pfam.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# sctc.s_peakedT3.pfam <- 
#   pfam %>% filter(query_name %in% sctc.s_peakedT3$orig) %>%
#   group_by(description) %>%
#   count(sort = TRUE, name = "count")

sctc.s_peakedT3.pfam <- pfam %>% 
  filter(query_name %in% sctc.s_peakedT3$orig) %>%
  select(query_name, domain_name, sequence_evalue, description)

str(sctc.s_peakedT3.pfam)
# tibble [62 × 4] (S3: tbl_df/tbl/data.frame)

write.table(sctc.s_peakedT3.pfam, 
            "results/SCTC.S_PeakedT3.pfam.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


# sctc.r_peakedT2.pfam <- 
#   pfam %>% filter(query_name %in% sctc.r_peakedT2$orig) %>%
#   group_by(description) %>%
#   count(sort = TRUE, name = "count")

sctc.r_peakedT2.pfam <- pfam %>% 
  filter(query_name %in% sctc.r_peakedT2$orig) %>%
  select(query_name, domain_name, sequence_evalue, description)

str(sctc.r_peakedT2.pfam)
# tibble [399 × 4] (S3: tbl_df/tbl/data.frame)

write.table(sctc.r_peakedT2.pfam, 
            "results/SCTC.R_PeakedT2.pfam.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)


# sctc.r_peakedT3.pfam <- 
#   pfam %>% filter(query_name %in% sctc.r_peakedT3$orig) %>%
#   group_by(description) %>%
#   count(sort = TRUE, name = "count")

sctc.r_peakedT3.pfam <- pfam %>% 
  filter(query_name %in% sctc.r_peakedT3$orig) %>%
  select(query_name, domain_name, sequence_evalue, description)

str(sctc.r_peakedT3.pfam)
# tibble [536 × 4] (S3: tbl_df/tbl/data.frame)

write.table(sctc.r_peakedT3.pfam, 
            "results/SCTC.R_PeakedT3.pfam.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
