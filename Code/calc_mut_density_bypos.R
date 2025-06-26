library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

setwd("/data/gtl/avaquer/tRNA_mut_density/")

#Loading the files

files = list.files("/data/gds/franklin/cancer3_backup/mmurillo/results/mut_rates_per_tRNA_control_trint_PERPOS/", pattern = "*", full.names = T)
length(files)

#Removing the files that we are not going to use
files = files[!grepl("10.csv", files)]
files = files[!grepl("NA.csv", files)]
length(files)
#Loading metadata

metadata = read_tsv("/data/gtl/avaquer/tRNA_mut_density/data/metadata_final.tsv")
metadata <- metadata %>%
  mutate(tissue = ifelse(tissue == "" | is.na(tissue) , "unknown", tissue)) %>%  # Replace NA or empty with "unknown"
  group_by(sample_id) %>%
  summarise(tissue = first(tissue))  # Collapse tissue by commas

# Creating a data frame to iterate trough the different files

df_files = data.frame(file_name = files)
head(df_files)
df_files$tRNA = sub(".*/", "", df_files$file_name)
df_files$tRNA = sub(".csv", "", df_files$tRNA)


# Iterating over the files to obtain the mutational density per position and per patient
df_counts = dplyr::bind_rows(lapply(df_files$file_name, function(f){
  d = read_csv(f)
  d = merge(d, metadata, by= "sample_id", all.x = TRUE)
  df_counts = d
  df_counts = df_counts %>% group_by(group, context) %>%
    summarise(sum_count_muts = sum(count_muts),sum_nt_risk = sum(count)) %>%
    mutate(tissue = ifelse(tissue == "" | is.na(tissue) , "unknown", tissue))
  df_counts$density = as.numeric(df_counts$sum_count_muts/df_counts$sum_nt_risk)
  df_counts_red = df_counts %>% group_by(group) %>% summarise(sum_density = sum(density, na.rm = TRUE)) %>% mutate(tRNA = sub("_.+", "", group),           # Extract part before the first underscore
                                                                                                                       pos = sub(".+__","", group)                   # Extract number after the last double underscore
                                                                                                                     )
  df_counts_red = df_counts_red[,c("tRNA", "pos", "sum_density")]
  head(df_counts_red)
  write.csv(df_counts_red, paste0("/data/gtl/avaquer/tRNA_mut_density/results_mut_density_bypos/sum_count_density_with_trint_by_pos",basename(f)), row.names = F)
}))





