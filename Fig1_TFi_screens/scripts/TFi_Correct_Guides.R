#Corrects duplicated guides from the TFi Library
#Some guides were designed twice, if they targeted nearby promoters associated with different genes
#This script adds back the lost assignments to the count tables
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
input_dir <- args[2] 
output_dir <- args[3]

setwd(output_dir)

lib_file <- read.delim(paste0(input_dir, "TFi_library_file.txt"),
                       col.names = c("sgRNA", "seq", "Gene"), header = FALSE)

counts_raw <- read.delim(paste0(data_dir, "TFi.count.txt"))
counts_norm <- read.delim(paste0(data_dir, "TFi.count_normalized.txt"))

counts_seq_raw <- counts_raw %>% 
  left_join(lib_file)

counts_seq_norm <- counts_norm %>% 
  left_join(lib_file)

#Checks for guides that are included in the library file, but not in the count tables 
#(MAGeCK removes one of the duplicated guides)
miss_guides <- lib_file %>% 
  filter(!sgRNA %in%  counts_raw$sgRNA) %>% 
  dplyr::rename(miss_sgRNA = sgRNA, miss_Gene = Gene) 

miss_guides_raw <- inner_join(miss_guides, counts_seq_raw) %>% 
  select(-seq, -sgRNA, -Gene) %>% 
  dplyr::rename(sgRNA = miss_sgRNA, Gene = miss_Gene)

miss_guides_norm <- inner_join(miss_guides, counts_seq_norm) %>% 
  select(-seq, -sgRNA, -Gene) %>% 
  dplyr::rename(sgRNA = miss_sgRNA, Gene = miss_Gene)

#Adds the lost assignments back to the count tables
total_counts_raw <- rbind(counts_raw, miss_guides_raw)
total_counts_norm <- rbind(counts_norm, miss_guides_norm)

#Writes the corrected count tables to the output
write_delim(total_counts_raw,
            paste0(output_dir, "TFi.count_corr.txt"), delim = "\t")

write_delim(total_counts_norm,
            paste0(output_dir, "TFi.count_corr_normalized.txt"), delim = "\t")
