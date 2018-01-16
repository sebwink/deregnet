#!/opt/anaconda/3.5/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

source("deseq2.R")

# TODO: Use optparse, do not forget to rewrite deseq2.py afterwards

count_matrix <- load_count_matrix(args[1])
coldata <- load_coldata(args[2])
design <- "~condition"
num_cores <- 4

cat("\n########## Normalizing counts with DESeq2 ##########\n\n")
normalized_counts <- get_normalized_counts(count_matrix, coldata, design, num_cores)
cat("\nDone with normalization. Writing results.\n")
write.csv(normalized_counts, "normalized_counts.csv")
