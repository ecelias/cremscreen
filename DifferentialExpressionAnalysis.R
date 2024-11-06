library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readr)
library(DESeq2)

merged_data_csv <- "merged_data_with_clusters.csv"
merged_counts <- read.csv(merged_data_csv)
names(merged_counts)[1] <- "CellID"

rownames(merged_counts) <- merged_counts$CellID
merged_counts$CellID <- NULL

diff_exp <- function(data){
  cluster_labels <- as.factor(data$Cluster)
  expression_data <- data[, -ncol(data)]
  expression_data <- expression_data[rowSums(expression_data) > 0, ]  # Remove Cluster column
  
  # Create a DESeqDataSet
  col_data <- data.frame(cluster = cluster_labels)
  dds <- DESeqDataSetFromMatrix(countData = t(expression_data), 
                                colData = col_data, 
                                design = ~ cluster)
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  results <- results(dds)
  
  # Order results by p-value and log2 fold change
  results <- results[order(results$padj, -abs(results$log2FoldChange)), ]
  
  # Get top 10 differentially expressed genes
  top_genes <- head(results, 10)
  print(top_genes)
  
}

diff_exp(merged_counts)


