library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(readr)

file1 <- "merged_atac_rna_counts.csv"
merged_counts <- read_csv(file1)