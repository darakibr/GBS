
### R script for hierbaps ###
library(rhierbaps)
library(ggtree)
library(phytools)

snp.matrix <- load_fasta("/home/user/Documents/GBS_typeIV/FOLDER/vcf-filtered/FOLDER_exclude.fa")

hb.results <- hierBAPS(snp.matrix, max.depth=2, n.pops=30, quiet=TRUE)

write.csv(hb.results$partition.df, file="/home/user/Documents/GBS_typeIV/FOLDER/FOLDER_hierbaps.csv", col.names=TRUE, row.names=FALSE)
  
