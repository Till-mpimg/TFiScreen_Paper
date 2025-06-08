#Filters FIMO motif enrichment for interesting occurences
library(tidyverse)



#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Reads RE and motif info
re_total <- read.delim("./input_files/re_total.bed",
                       col.names = c("chr", "start", "end", "RE", "trash1", "trash2", "trash3", "trash4", "trash5"),
                       header = FALSE) %>% 
  transmute(sequence_name = paste0(chr, ":", start, "-", end), RE)

motifs <- read.delim("./output_files/TFiMini_motifs.txt", sep = "\t",
                     header = FALSE, col.names = c("motif_id", "cluster_id", "domain", "tf", "family", "source", "pmid")) %>% 
  select(cluster_id, tf) %>% 
  unique()

motifs_conc <- motifs %>% 
  group_by(cluster_id) %>% 
  summarise(tfs_conc = paste(unique(tf), collapse = ';')) %>% 
  mutate(name = paste0(cluster_id, " (", tfs_conc, ")"))

#Reads FIMO Info
fimo <- read.delim("./output_files/fimo.tsv") %>% 
  na.omit() %>% 
  mutate(cluster_id = str_extract(motif_id, "AC[0-9]*")) %>% 
  left_join(re_total) %>% 
  inner_join(motifs_conc, multiple = "all")

fimo_count <- fimo %>% 
  select(name, RE, score) %>% 
  group_by(name, RE) %>% 
  mutate(RE = str_remove(RE, "_"), name = str_replace(name, "ZNF", "ZFP")) %>% 
  summarize(n = n(), max_score = max(score))

RE_vec <- c("RE12", "RE46", "RE47", "RE49", "RE50", "RE51", "RE52", "RE53", "RE55", 
            "RE57L", "RE57M", "RE57R", "RE58", "RE59", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")


#Writes info to file
write_delim(fimo_count, "./output_files/motifs_re_clean.txt", delim = "\t")
