#!/usr/bin/env Rscript
library(plyr)
args = commandArgs(trailingOnly=TRUE)
 genes <-data.frame(read.table(args[1]))
names <-data.frame(read.table(args[2]))
names(genes)<- c("gene", "short.name")
head(genes)
names(names)<- c("gene", "short.name")
head(names)
merged_and_ordered <- join(x=genes, y=names, by="gene", type = "left")

write.table(merged_and_ordered, file=args[3])


