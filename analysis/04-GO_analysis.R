library(dplyr)
library(purrr)
library(tibble)
library(tidyr)


sigDE <- read.table("results/StoneCellBiosyn_pooledRun.sigDE.25oct.txt", 
                       header = TRUE, stringsAsFactors = FALSE)
str(sigDE)
# 'data.frame':	6565 obs. of  5 variables:


cType_DE <- subset(sigDE, sigDE$logFC > 0 & 
                     (sigDE$focus == "H898_cType" | sigDE$focus == "Q903_cType"))
cType_DE <- cType_DE %>% select("contig", "focus")
str(cType_DE)
# 'data.frame':	1021 obs. of  2 variables:

interx <- cType_DE[(duplicated(cType_DE$contig)), "contig"]
length(interx)
# [1] 70

h898_excl <- cType_DE %>% filter(!contig %in% interx & focus == "H898_cType") %>% 
  mutate(set = "H898_excl") %>% select(contig, set) %>% distinct()
intersection <- cType_DE %>% filter(contig %in% interx) %>% mutate(set = "intersection") %>% 
  select(contig, set) %>% distinct()
q903_excl <- cType_DE %>% filter(!contig %in% interx & focus == "Q903_cType") %>% 
  mutate(set = "Q903_excl") %>% select(contig, set) %>% distinct()


cType_DE <- rbind((rbind(h898_excl, intersection)), q903_excl)
dim(cType_DE)
# [1] 951   2

table(cType_DE$set)
#    H898_excl intersection    Q903_excl 
#          738           70          143 


# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' 
# option on comamnd line BLAST).  We ran with 10 max target hits returned 
# per query sequences.
TAIR_Annots <- read.delim("data/cType_DE.blastpTAIR10.txt", 
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
TAIR_Annots <- TAIR_Annots %>% select(qseqid, sseqid, evalue)
str(TAIR_Annots)
# 'data.frame':	8339 obs. of  3 variables:


# There are as much as 10 subject sequence per ID.  We are only using the top
# hit for the purpose of pathway annotation
getTopHit <- function(annots) {
  topHit <- annots %>% group_by(qseqid) %>% nest()
  topHit$data <- map(topHit$data, function(x) {x[1,]})
  topHit <- unnest(topHit)
}

top_TAIR_Annots <- getTopHit(TAIR_Annots)
str(top_TAIR_Annots)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	1039 obs. of  3 variables:


# Merge DE results with top TAIR annotation
cType_DE_TAIR_Annots <- merge(cType_DE, top_TAIR_Annots, by.x = "contig", by.y = "qseqid", all.x = TRUE)
# 'data.frame':	6565 obs. of  6 variables:


### Read Ath (TAIR10) gene ontology
geneOnt <- read.delim("data/ATH_GO_GOSLIM.txt", header = FALSE, 
                      sep = "\t", stringsAsFactors = FALSE)
str(geneOnt)
# 'data.frame':	336022 obs. of  15 variables:
 
# From TAIR readme file
# 1. locus name: standard AGI convention name
# 2. TAIR accession:the unique identifier for an object in the TAIR database- 
#   the object type is the prefix, followed by a unique accession number(e.g. gene:12345).  
# 3. object name : the name of the object (gene, protein, locus) being annotated.
# 4. relationship type: the relationship between the annotated object and the GO term
# 5. GO term: the actual string of letters corresponding to the GO ID
# 6. GO ID: the unique identifier for a GO term.  
# 7. TAIR Keyword ID: the unique identifier for a keyword in the TAIR database.
# 8. Aspect: F=molecular function, C=cellular component, P=biological 13process. 
# 9. GOslim term: high level GO term helps in functional categorization.
# 10. Evidence code: three letter code for evidence types 
# 11. Evidence description: the analysis that was done to support the annotation
# 12. Evidence with: supporting evidence for IGI, IPI, IC, IEA and ISS annotations
# 13. Reference: Either a TAIR accession for a reference
# 14. Annotator: TAIR, TIGR or a TAIR community member
# 15. Date annotated: date the annotation was made.

geneOnt <- geneOnt %>% select(3,5,8)
colnames(geneOnt) <- c("locus", "GO_term", "Aspect")


geneOnt_c <- subset(geneOnt, geneOnt$Aspect == "C")
geneOnt_c <- unique(geneOnt_c)
dim(geneOnt_c)
# [1] 62298     3

geneOnt_c <- geneOnt_c %>% select(-Aspect) %>% group_by(locus) %>% nest()
geneOnt_c$data <- map(geneOnt_c$data, function(x){
  paste0(unlist(x), collapse = " ; ")
})
geneOnt_c <- geneOnt_c %>% unnest()
str(geneOnt_c)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	42812 obs. of  2 variables:


geneOnt_f <- subset(geneOnt, geneOnt$Aspect == "F")
geneOnt_f <- unique(geneOnt_f)
dim(geneOnt_f)
# [1] 66764     3

geneOnt_f <- geneOnt_f %>% select(-Aspect) %>% group_by(locus) %>% nest()
geneOnt_f$data <- map(geneOnt_f$data, function(x){
  paste0(unlist(x), collapse = " ; ")
})
geneOnt_f <- geneOnt_f %>% unnest()
str(geneOnt_f)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	41580 obs. of  2 variables:


geneOnt_p <- subset(geneOnt, geneOnt$Aspect == "P")
geneOnt_p <- unique(geneOnt_p)
dim(geneOnt_p)
# [1] 89350     3

geneOnt_p <- geneOnt_p %>% select(-Aspect) %>% group_by(locus) %>% nest()
geneOnt_p$data <- map(geneOnt_p$data, function(x){
  paste0(unlist(x), collapse = " ; ")
  })
geneOnt_p <- geneOnt_p %>% unnest()
str(geneOnt_p)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	38263 obs. of  2 variables:


results <- merge(cType_DE_TAIR_Annots, geneOnt_c, by.x = "sseqid", by.y = "locus", all.x = TRUE)
results <- merge(results, geneOnt_f, by.x = "sseqid", by.y = "locus", all.x = TRUE)
results <- merge(results, geneOnt_p, by.x = "sseqid", by.y = "locus", all.x = TRUE)
str(results)
# 'data.frame':	6565 obs. of  9 variables:

colnames(results) <- c("TAIR_id", "contig", "set", "evalue", 
                       "GO_Cell.Comp", "GO_Mol.Func", "GO_Biol.Proc")
results <- 
  results %>% select("contig", "TAIR_id", "set", "evalue", 
                     "GO_Cell.Comp", "GO_Mol.Func", "GO_Biol.Proc")

write.table(results, "results/StoneCellBiosyn.cType_DE.GO_term.txt", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
