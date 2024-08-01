#Script to cluster TFs according to their expression profiles
library(stats)
library(tidyverse)
library(EnvStats)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


#Reads in data from TX timecourse and TFi Screen
cpm_tx <- read.delim("./output_files/CPM_RNA_timecourse.txt")
tf_list <- read.delim("./input_files/Mus_musculus_TF.txt")$Symbol #From TF database
lib_file <- read.delim("./input_files/TFi_isoforms_clean.txt")

TFi_TFs <- lib_file %>% 
  filter(gene %in% tf_list) %>%
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

#Creates table with expression and only focuses on TFs and Xist
cpm_long <- cpm_tx %>%  
  filter(gene %in% TFi_TFs) %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("sex", "rep", "day"))
  

#make this example reproducible
set.seed(11)

max_day <- cpm_long %>% 
  group_by(gene, gene_id, day) %>% 
  summarize(cpm = mean(cpm)) %>% 
  ungroup() %>% 
  group_by(gene, gene_id) %>% 
  mutate(max_day = day[which.max(cpm)]) %>% 
  select(gene, gene_id, max_day)

scale_cpm <- cpm_long %>% 
  filter(gene != "Xist") %>% 
  group_by(gene_id) %>% 
  mutate(zscore = as.numeric(scale(cpm)))

pre_kmeans <- scale_cpm %>% 
  left_join(max_day) %>% 
  group_by(sex, gene, gene_id, day, max_day) %>% 
  summarize(zscore = mean(zscore))

#Cluster Formative TFs
form_wide <- pre_kmeans %>% 
  ungroup() %>% 
  filter(max_day %in% c("10h", "16h", "24h", "30h")) %>% 
  unite(sample, sex, day) %>% 
  select(-gene, -max_day) %>% 
  pivot_wider(names_from = sample, values_from = zscore) %>%  
  column_to_rownames("gene_id")

kmeans_form <- kmeans(form_wide , 2)

#Cluster Primed TFs
prime_wide <- pre_kmeans %>% 
  ungroup() %>% 
  filter(max_day %in% c("36h", "48h", "56h", "72h", "96h")) %>% 
  unite(sample, sex, day) %>% 
  select(-gene, -max_day) %>% 
  pivot_wider(names_from = sample, values_from = zscore) %>%  
  column_to_rownames("gene_id")


kmeans_prime <- kmeans(prime_wide, 2)


#Combine clusters
cluster_form <- data.frame(cluster = paste0("form", kmeans_form$cluster), gene_id = names(kmeans_form$cluster))
cluster_prime <- data.frame(cluster = paste0("prime", kmeans_prime$cluster), gene_id = names(kmeans_prime$cluster))
cluster_naive <- pre_kmeans %>% 
  ungroup() %>% 
  filter(max_day == "0h") %>% 
  select(gene_id) %>%
  unique() %>% 
  mutate(cluster = "naive")

clusters <- rbind(cluster_naive, cluster_form, cluster_prime)


#Give out clusters
write_delim(clusters, "./output_files/Timecourse_TF_clusters.txt", delim = "\t")
