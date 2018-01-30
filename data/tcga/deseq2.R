library(BiocParallel)
library(DESeq2)

load_count_matrix <- function(path2countMatrix) {
  # ---
  count_matrix <- read.table(gzfile(path2countMatrix), header = TRUE, check.names = FALSE)
  rownames(count_matrix) <- count_matrix[,"id"]
  count_matrix <- subset(count_matrix, select=-c(id))
  count_matrix <- count_matrix[, sort(colnames(count_matrix))] 
  return(count_matrix)
}

load_coldata <- function(path2coldata) {
  # ---
  coldata <- read.table(gzfile(path2coldata), header = TRUE, check.names = FALSE)
  rownames(coldata) <- coldata[,"file_id"]
  coldata <- coldata[sort(rownames(coldata)),]
  coldata <- subset(coldata, select=-c(file_id))
  return(coldata)
}

get_normalized_counts <- function(count_matrix,
								  coldata,
								  design,
								  num_cores = 4) {
  register(MulticoreParam(num_cores))
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = coldata,
                                design = as.formula(design))
  dds$conditions <- factor(dds$condition,levels = c("normal", "tumor")) # ?
  dds=dds[rowSums(counts(dds))>1,]
  dds <- estimateSizeFactors(dds)
  return(counts(dds, normalized=TRUE))
}
