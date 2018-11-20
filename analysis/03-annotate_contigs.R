library(dplyr)

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


annot <- read.table("data/cType_DE.blastpNR.txt", header = TRUE, sep = "\t", quote = "",
                    stringsAsFactors = FALSE)
annot <- annot[ ,c("qseqid", "salltitles")]
str(annot)
# 'data.frame':	12940 obs. of  2 variables:


concatAnnot <- function(x) {
  title <- x$salltitles
  # title <- gsub('<>', " ", annot)
  return (paste0(title, collapse = " ; "))
}


annotsCollapsed <-
  ddply(annot, ~ qseqid, concatAnnot)
colnames(annotsCollapsed) <- c("contig", "annot")
str(annotsCollapsed)
# 'data.frame':	1321 obs. of  2 variables:


merged <- merge(cType_DE, annotsCollapsed, by.x = "contig", by.y = "contig", all.x = TRUE)
str(merged)
# 'data.frame':	951 obs. of  3 variables:


write.table(merged, file = "results/StoneCellBiosyn.cType_DE.annot.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
