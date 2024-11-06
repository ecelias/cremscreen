# Load libraries
library(DESeq2)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readr)

# Load data
merged_data <- read.csv("merged_data_with_clusters.csv")

# Separate clusters
cluster_data <- merged_data$Cluster
counts_data <- merged_data
counts_data$Cluster <- NULL

counts_data_t <- transpose(counts_data)
rownames(counts_data_t) <- colnames(counts_data)
colnames(counts_data_t) <- rownames(counts_data)

colnames(counts_data_t) <- counts_data_t[1,]
counts_data_t <- counts_data_t[-1,]

rownames(counts_data) <- counts_data[,1]
counts_data <- counts_data[,-1]

data_matrix <- as.matrix(counts_data_t)
clusters <- data.frame(cluster = as.factor(cluster_data))
clusters_t <- transpose(clusters)


# Define the function to perform differential expression analysis
differential_expression_analysis <- function(expression, meta) {
  
  data_matrix <- data.matrix(expression)
  clusters <- as.data.frame(meta)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                                colData = clusters,
                                design = ~ cluster)
  
  
  # Perform DE analysis, with workaround for dispersion error
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)  # Use gene-wise estimates for dispersion
  dispersions(dds) <- mcols(dds)$dispGeneEst  # Set dispersions
  dds <- nbinomWaldTest(dds)  # Continue with testing
  
  res <- results(dds)
  
  # Sort results by adjusted p-value and select top 10 genes
  res <- res[order(res$padj), ]
  top_genes <- head(res, 10)
  
  # Print top genes
  print(top_genes)
  
  # Plot volcano plot
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  title = 'Volcano plot of differential expression',
                  pCutoff = 0.05,
                  FCcutoff = 1)
}

# Run the function
differential_expression_analysis(expression=counts_data, meta=clusters_t)



