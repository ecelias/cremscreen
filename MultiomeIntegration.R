library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readr)

plot_multiome <- function(multiome_sum) {
  multiome_comparison <- ggplot(multiome_sum, aes(x=merged_frag_sums, y=merged_gene_sums)) + geom_point()
  ggsave(multiome_comparison, device = "png")
}

multiome_analysis <- function(filename1, filename2){
  # read in the csv files
  total_counts <- read_csv(filename1)
  fragment_counts <- read.csv(filename2)
  
  # transpose the data frames (using data.table)
  frag_t <- transpose(fragment_counts)
  total_t <- transpose(total_counts)
  
  # update row and column names back to original 
  rownames(frag_t) <- colnames(fragment_counts)
  colnames(frag_t) <- rownames(fragment_counts)
  rownames(total_t) <- colnames(total_counts)
  colnames(total_t) <- rownames(total_counts)
  
  colnames(frag_t) <- frag_t[1,]
  frag_t <- frag_t[-1,]
  colnames(total_t) <- total_t[1,]
  total_t <- total_t[-1,]
  
  # create a copy of original data frames but
  # add rows as a new column
  copy_frag_t <- frag_t   
  copy_frag_t <- tibble::rownames_to_column(copy_frag_t, "CellID") 
  
  copy_total_t <- total_t   
  copy_total_t <- tibble::rownames_to_column(copy_total_t, "CellID") 
  
  # merge with inner_join (dyplr) to remove any cell IDs not found in both datasets
  merged <- inner_join(copy_frag_t, copy_total_t, by="CellID")
  
  # updated the merged dataframes row names as the "CellID" column
  rownames(merged) <- merged$CellID
  merged$CellID <- NULL
  
  # convert the data frame to a numeric matrix
  merged_matrix <- data.matrix(merged)
  
  # get a vector of columns for either chromosome regions or genes
  chr_columns <- grep("^chr", colnames(merged_matrix))
  gene_columns <- grep("^Gene", colnames(merged_matrix))
  
  # find the sum of the regions and genes for each cell 
  merged_frag_sums <- rowSums(merged_matrix[,chr_columns], na.rm=TRUE)
  merged_gene_sums <- rowSums(merged_matrix[,gene_columns], na.rm=TRUE)
  
  sums_df <- data.frame(
    fragments = c(merged_frag_sums),
    genes = c(merged_gene_sums)
  )
  
  plot_multiome(sums_df)
}

# sample usage
filename1 <- 'scATAC_counts.csv'
filename2 <- 'scATAC_fragment_counts.csv'
multiome_analysis(filename1, filename2)


