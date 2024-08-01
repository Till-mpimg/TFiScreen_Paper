#Script to perform DESeq2 on TX timecourse data
library(tidyverse)
library(DESeq2)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Calculate DESEQ2
counts <- read.delim("./output_files/Counts_RNA_timecourse.txt") 

time_vec <- c("0h", "10h", "16h", "24h", "30h", "36h", "48h", "56h", "72h", "96h")

deseq_fun <- function(timepoint) {
  fil_counts <- counts %>% 
    pivot_longer(-c(gene, gene_id), names_to = "sample", values_to = "count") %>% 
    filter(str_detect(sample, paste0("_", timepoint))) %>% 
    arrange(sample)
  
  fil_wide <- fil_counts %>% 
    pivot_wider(names_from = sample, values_from = count) %>%
    select(-gene) %>% 
    column_to_rownames("gene_id")
  
  line <- str_extract(unique(fil_counts$sample), "X[XO]")
  time_vec <- rep(timepoint, length(line))
  coldata <- data.frame(line, time_vec)
  rownames(coldata) <- colnames(fil_wide)
  
  dds <- DESeqDataSetFromMatrix(countData = fil_wide,
                                colData = coldata,
                                design = ~ line)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  res <- data.frame(results(DESeq(dds))) %>%
    rownames_to_column("gene_id") %>% 
    mutate(timepoint = timepoint)
  return(res)
}

#Puts all DESEQ comparisons into a single dataframe and exports to a text file
deseq_total <- bind_rows(lapply(time_vec, deseq_fun))
deseq_df <- deseq_total %>% 
  left_join(counts[,c(1:2)])

write_delim(deseq_df, "./output_files/DEseq2_total.txt", delim = "\t")
