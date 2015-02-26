library(edgeR)

source("analysis/helper02-DE_workflow.R")
source("analysis/helper04-proc_results.R")

run <- "julyRun"
julyFit <- DE_workflow(run, TRUE)

# N.B.  LogFC cut off threashold set at 2
summary(decideTests(julyFit, method = "global", adjust.method = "fdr", 
                    p.value = 0.01, lfc = 2))
# For July Assembly:
#    CP_gType DSC_gType Q903_cType H898_cType
# -1      837       871         77        119
# 0     34195     33894      35651      35364
# 1       772      1039         76        321

julyResults <- processResults(julyFit)

summary(julyResults$focus)
# For July Assembly:
#   CP_gType  DSC_gType Q903_cType H898_cType 
#      35804      35804      35804      35804 

julyResults.wide <- reshape(julyResults, direction = "wide", 
                            idvar = "contig", timevar = "focus")

### Write results files
write.table(
  julyResults, "results/StoneCellBiosyn_julyRun_DE_results.long.txt", 
  quote = FALSE, sep = "\t", row.names = FALSE)

write.table(
  julyResults.wide, "results/StoneCellBiosyn_julyRun_DE_results.wide.txt",
  quote = FALSE,  sep = "\t", row.names = FALSE)

#####

run <- "novRun"
novFit <- DE_workflow(run, TRUE)

summary(decideTests(novFit, method = "global", adjust.method = "fdr", 
                    p.value = 0.01, lfc = 2))
# For Nov Assembly:
#    CP_gType DSC_gType Q903_cType H898_cType
# -1     1127      1068        163        142
# 0     32349     32277      34156      33903
# 1       931      1062         88        362

novResults <- processResults(novFit)

summary(novResults$focus)
# For Nov Assembly:
#   CP_gType  DSC_gType Q903_cType H898_cType 
#      34407      34407      34407      34407

novResults.wide <- reshape(novResults, direction = "wide", 
                           idvar = "contig", timevar = "focus")

### Write results files
write.table(
  novResults, "results/StoneCellBiosyn_novRun_DE_results.long.txt", 
  quote = FALSE, sep = "\t", row.names = FALSE)

write.table(
  novResults.wide, "results/StoneCellBiosyn_novRun_DE_results.wide.txt",
  quote = FALSE,  sep = "\t", row.names = FALSE)


### Write out all DE contig ids with abs(logFC) >= 2 & adj. p-value <= 0.01
write.table(
  subset(
    julyResults, julyResults$global.adj.P <= 0.01 & abs(julyResults$logFC) >= 2)[, "contig"],
  "results/StoneCellBiosyn_julyRun_DE_lfc2.id.txt.new", quote = FALSE,
  sep = "\t", col.names = FALSE, row.names = FALSE)

write.table(
  subset(novResults, novResults$global.adj.P <= 0.01 & abs(novResults$logFC) >= 2)[, "contig"],
  "results/StoneCellBiosyn_novRun_DE_lfc2.id.txt.new", quote = FALSE,
  sep = "\t", col.names = FALSE, row.names = FALSE)
