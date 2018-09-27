library(plyr)

sigDE <- read.table("results/StoneCellBiosyn_pooledRun_sigDE.18aug.txt", header = TRUE)
str(sigDE)
# 'data.frame':	17913 obs.


annot <- read.table("results/StoneCellBiosyn_pooledRun_sigDE.18aug.fasta.blastxNr.txt", 
                    header = TRUE, sep = "\t", quote = "")
annot <- annot[,c(1,25)]
str(annot)
# 'data.frame':	50578 obs. of  5 variables:


concatAnnot <- function(x) {
  title <- x$salltitles
  # title <- gsub('<>', " ", annot)
  return (paste0(title, collapse = " ; "))
}


annotsCollapsed <-
  ddply(annot, ~ qseqid, concatAnnot)
colnames(annotsCollapsed) <- c("contig", "annot")


merged <- merge(sigDE, annotsCollapsed, by.x = "contig", by.y = "contig", all.x = TRUE)


write.table(merged, file = "results/StoneCellBiosyn_pooledRun_sigDE.18aug.annot.txt",
            row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
