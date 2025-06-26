library(readr)
library(dplyr)
library(MutationalPatterns)
library(tibble)
library(patchwork)
library(stringr)
library(viridis)
library(ggplot2)
# Loading the data and keeping tRNA and protein-coding data separately
mut_data <- read.csv("/Users/ainav/OneDrive/Escritorio/Prácticas IRB/MutSignatures/df_ms96.csv")
meta <-read_tsv("/Users/ainav/OneDrive/Escritorio/Prácticas IRB/MutSignatures/comb_metadata_final_6datasets_updatedNov2021_corrected.txt", show_col_types = FALSE)
mut_by_tiss <- read.csv("/Users/ainav/OneDrive/Escritorio/Prácticas IRB/MutSignatures/mutsig_by_tiss.csv", sep =";")
meta$tissue[meta$tissue=="colorectum" & meta$MSI_status == "MSI"] = "colorectum_MSI"
meta$tissue[meta$tissue=="colorectum" & meta$MSI_status == "MSS"] = "colorectum_MSS"
meta$tissue[meta$tissue=="colorectum" & meta$MSI_status %in% c("ERROR", "HYPER")] = "colorectum_HYPER"
meta$tissue[meta$tissue=="uterus" & meta$MSI_status == "MSI"] = "uterus_MSI"
meta$tissue[meta$tissue=="uterus" & meta$MSI_status == "MSS"] = "uterus_MSS"
meta$tissue[meta$tissue=="uterus" & meta$MSI_status %in% c("ERROR", "HYPER")] = "uterus_HYPER"
meta <- meta[,c("tissue")] %>% group_by(tissue) %>% summarise(nsamples = n())
meta$tissue <- str_replace(meta$tissue, "^[^_]+", function(x) {
  str_c(str_to_upper(str_sub(x, 1, 1)), str_sub(x, 2))
})
mut_data$tissue <- str_replace(mut_data$tissue, "^[^_]+", function(x) {
  str_c(str_to_upper(str_sub(x, 1, 1)), str_sub(x, 2))
})
mut_data$tissue <- gsub("_", " ", mut_data$tissue)
meta$tissue <- gsub("_", " ", meta$tissue)
meta <- meta[meta$tissue %in% mut_by_tiss$Tissue,]
mut_data_tRNA <- merge(mut_data[mut_data$i == 0 & mut_data$gene_type == "tRNA", ], meta, by = 'tissue')  %>%
  mutate(mut_type = gsub("\\(", "[", mut_type),
         mut_type = gsub("\\)", "]", mut_type)) %>%
  group_by(mut_type, tissue) %>%
  summarise(mut_count = sum(sum_count_muts))
mut_data_prot <- merge(mut_data[mut_data$i == 0 & mut_data$gene_type == "proteincoding", ], meta, by = 'tissue') %>%
  mutate(mut_type = gsub("\\(", "[", mut_type),
         mut_type = gsub("\\)", "]", mut_type)) %>%
  group_by(mut_type, tissue) %>%
  summarise(mut_count = sum(sum_count_muts))
# Selecting unique tissues to use as rows of the mut matrix table
tissues <- unique(mut_data_tRNA$tissue)
#Getting all possible contexts from COSMIC
# Define bases and substitutions
bases <- c("A", "C", "G", "T")
substitutions <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
# Generate all trinucleotide combinations
contexts <- c()
for (sub in substitutions) {
  for (first in bases) {
    for (third in bases) {
      contexts <- c(contexts, paste0(first, "[", sub, "]", third))
    }
  }
}
cosmic_sig <- get_known_signatures(muttype = "snv", genome = 'GRCh37')
mut_matrix_tRNA <-  data.frame(matrix(0, nrow = 96, ncol = length(tissues )))
colnames(mut_matrix_tRNA) <- tissues
rownames(mut_matrix_tRNA) <- contexts
mut_matrix_prot <-  data.frame(matrix(0, nrow = 96, ncol = length(tissues)))
colnames(mut_matrix_prot) <- tissues
rownames(mut_matrix_prot) <- contexts
# Fill in mut_matrix for tRNA genes
for (i in 1:nrow(mut_data_tRNA)) {
  row <- mut_data_tRNA$mut_type[i]
  col <- mut_data_tRNA$tissue[i]
  count <- mut_data_tRNA$mut_count[i]
  # Check if row and column exist
  if (row %in% rownames(mut_matrix_tRNA) && col %in% colnames(mut_matrix_tRNA)) {
    mut_matrix_tRNA[row, col] <- mut_matrix_tRNA[row, col] + count
  } else {
    warning(paste("Skipping unknown context or tissue:", row, col))
  }
}
# Fill in mut_matrix for protein-coding
for (i in 1:nrow(mut_data_prot)) {
  row <- mut_data_prot$mut_type[i]
  col <- mut_data_prot$tissue[i]
  count <- mut_data_prot$mut_count[i]
  # Check if row and column exist
  if (row %in% rownames(mut_matrix_prot) && col %in% colnames(mut_matrix_prot)) {
    mut_matrix_prot[row, col] <- mut_matrix_prot[row, col] + count
  } else {
    warning(paste("Skipping unknown context or tissue:", row, col))
  }
}

#Creating empty dataframes for tRNA and protein coding genes
signatures_fitted_tRNA <- data.frame(tissues = character(0))
# Add numeric columns for each signature (empty numeric vectors)
signatures_fitted_tRNA <- data.frame(matrix(ncol = length(colnames(cosmic_sig)), nrow = 0))
colnames(signatures_fitted_tRNA) <- colnames(cosmic_sig)
signatures_fitted_prot <- signatures_fitted_tRNA
# Obtain mutational signatures for tRNA
for(num in 1:ncol(mut_matrix_tRNA)){
  matrix_subset <- as.data.frame(mut_matrix_tRNA[, num, drop = FALSE])
  found_signatures <- mut_by_tiss[mut_by_tiss$Tissue == colnames(matrix_subset), ]
  found_signatures <- names(found_signatures)[which(found_signatures[1, ] == "Yes")]
  cosmic_subset <- cosmic_sig[,intersect(found_signatures, colnames(cosmic_sig))]
  temp_fit <- fit_to_signatures_bootstrapped(matrix_subset, cosmic_subset, method = 'strict', n_boots = 100)
  for (row in 1:nrow(temp_fit)) {
    row_name <- rownames(temp_fit)[row]
    zero_row <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(signatures_fitted_tRNA)))
    colnames(zero_row) <- colnames(signatures_fitted_tRNA)
    rownames(zero_row) <- row_name
    signatures_fitted_tRNA <- rbind(signatures_fitted_tRNA, zero_row)
    new_row_index <- nrow(signatures_fitted_tRNA)
    for (col in colnames(temp_fit)) {
      signatures_fitted_tRNA[new_row_index, col] <- temp_fit[row, col]
    }
  }
}

p <- plot_bootstrapped_contribution(signatures_fitted_tRNA,
                               mode = "relative",
                               plot_type = "dotplot")
p + scale_color_viridis_c(option = "plasma") +
  labs(y = "Tissue")


for(num in 1:ncol(mut_matrix_prot)){
  matrix_subset <- as.data.frame(mut_matrix_prot[, num, drop = FALSE])
  found_signatures <- mut_by_tiss[mut_by_tiss$Tissue == colnames(matrix_subset), ]
  found_signatures <- names(found_signatures)[which(found_signatures[1, ] == "Yes")]
  cosmic_subset <- cosmic_sig[,intersect(found_signatures, colnames(cosmic_sig))]
  temp_fit <- fit_to_signatures_bootstrapped(matrix_subset, cosmic_subset, method = 'strict', n_boots = 100)
  for (row in 1:nrow(temp_fit)) {
    row_name <- rownames(temp_fit)[row]
    zero_row <- as.data.frame(matrix(0, nrow = 1, ncol = ncol(signatures_fitted_prot)))
    colnames(zero_row) <- colnames(signatures_fitted_prot)
    rownames(zero_row) <- row_name
    signatures_fitted_prot <- rbind(signatures_fitted_prot, zero_row)
    new_row_index <- nrow(signatures_fitted_prot)
    for (col in colnames(temp_fit)) {
      signatures_fitted_prot[new_row_index, col] <- temp_fit[row, col]
    }
    if ("SBS5" %in% colnames(signatures_fitted_tRNA)) {
      signatures_fitted_prot[new_row_index, "SBS5"] <- signatures_fitted_prot[new_row_index, "SBS5"] * 0.2  
    }
    if ("SBS3" %in% colnames(signatures_fitted_tRNA)) {
      signatures_fitted_prot[new_row_index, "SBS3"] <- signatures_fitted_prot[new_row_index, "SBS3"] * 0.2  
    }
  }
}

plot_bootstrapped_contribution(signatures_fitted_prot,
                               mode = "relative",
                               plot_type = "dotplot")


