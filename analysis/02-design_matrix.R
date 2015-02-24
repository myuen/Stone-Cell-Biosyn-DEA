## write design matrix
expDes <- data.frame(sample = libNames,
                     gType = factor(c(rep("H898", 6), rep("Q903", 6)),
                                    levels = c("Q903", "H898")),
                     cType = factor(rep(c(rep("CP", 3), rep("DSC", 3)), 2),
                                    levels = c("CP", "DSC")),
                     bioRep = as.numeric(rep(c(1, 2, 3), 4)))
expDes$group <- with(expDes, interaction(gType, cType))
str(expDes)
expDes

write.table(expDes, "data/stone_cell_biosyn_exp_design.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
