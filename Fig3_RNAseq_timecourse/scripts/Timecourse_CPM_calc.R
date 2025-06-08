#Script calculating CPM from the timecourse RNA-seq data
library(tidyverse)
library(Rsubread)
library(Rgb)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Specifies path to BAM files
RNA_dir <-  paste0(wd, "data/final_bam/")
gencode_xert <- paste0(wd, "input_files/GENCODE_vM25_plus_Xert.gtf")

gene_names <- read.gtf(gencode_xert) %>%
  select(GeneID = gene_id, gene = gene_name) %>%
  unique()

#Count reads in final BAM files
setwd(RNA_dir)

temp = list.files(pattern="*.bam$", fill.names = TRUE)


feature_counts <- featureCounts(temp, annot.ext = gencode_xert,isGTFAnnotationFile = TRUE, isPairedEnd = FALSE,
                                GTF.featureType = "exon", strandSpecific = 1, allowMultiOverlap = TRUE)
counts <- data.frame(feature_counts$counts, feature_counts$annotation) %>%
  select(-Chr, -Start, -End, -Strand)
names(counts) <- gsub(x = names(counts), pattern = "RNA_", replacement = "")
names(counts) <- gsub(x = names(counts), pattern = "_sorted.bam", replacement = "")

#Calculates total counts and cpm
total_counts <- counts %>% 
  select(1:59) %>% 
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>% 
  group_by(sample) %>% 
  summarize(total_counts = sum(count))
  

CPM <- counts %>%
  select(-Length) %>%
  pivot_longer(-GeneID, names_to = "sample", values_to = "count") %>% 
  left_join(total_counts) %>% 
  mutate(CPM = (count * 1000000)/ total_counts) %>%
  select(-count, -total_counts) %>% 
  pivot_wider(names_from = sample, values_from = CPM) %>% 
  inner_join(gene_names) %>% 
  select(gene, gene_id = GeneID, everything())

setwd(wd)

count_out <- counts %>%  
  inner_join(gene_names) %>% 
  select(gene, gene_id = GeneID, everything()) %>% 
  select(-Length)

#Writes CPM and raw counts to output_files
write_delim(count_out,  "./output_files/Counts_RNA_timecourse.txt", 
            delim = "\t")
write_delim(CPM, "./output_files/CPM_RNA_timecourse.txt", 
            delim = "\t")

  
  
  
  
  


