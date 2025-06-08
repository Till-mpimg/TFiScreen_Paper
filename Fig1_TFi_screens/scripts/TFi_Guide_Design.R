#This script was used to produce guide sequences for TFi library cloning from Guidescan2
#Guides were produced with the GuideScan2 webtool for mm39 
library(tidyverse)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


#Creates dataframe to match mm10 with mm39 (Produced with LiftOver webtool from the UCSC Genome Browser)
guidescan_input_mm39 <- read.delim("./input_files/TFi_guidescan_500bp_mm39.bed",
                                   header = FALSE, col.names = c("chr", "mm39_start", "mm39_end", "id")) %>% 
  mutate(Region.name = paste0(chr, ":", mm39_start, "-", mm39_end))


#Reads guide locations from Guidescan in mm10 (File produced with LiftOver webtool from the UCSC Genome Browser)
guidescan_bed_mm10 <- read.delim("./intput_files/TFi_guidescan_locations_mm10.txt",
                                 col.names = "location", header = FALSE) %>% 
  mutate(chr = str_extract(location, "chr[X,0-9]*"), gRNA_start = str_extract(location, "(?<=:)[0-9]*"),
         gRNA_end = str_extract(location, "[0-9]*$"), num  = row_number())

#Reads the Guidescan results in mm39 and writes them to mm10
guidescan_results <- read.csv("./input_files/TFi_guidescan_results_mm39.csv") %>% 
  mutate(num = row_number()) %>% 
  left_join(guidescan_bed_mm10) %>% 
  left_join(guidescan_input_mm39) %>% 
  select(-Region.name) %>% 
  inner_join(guidescan_df)

guidescan_n <- guidescan_results %>% 
  group_by(id) %>% 
  summarize(n = n())

#Shows transcripts with less than 12 guides
missing_transcripts <- subset(guidescan_n, n < 12)$id

missing_transcripts_3 <- subset(guidescan_n, n < 4)$id

#Filters top 12 guides according to efficiency score
#Removes target promoters with less than 4 possible sgRNAs
guidescan_filtered <- guidescan_results %>%
  filter(!id %in% missing_transcripts_3) %>% 
  group_by(id) %>% 
  slice_max(Cutting.efficiency, n = 12)

#Creates sgRNA file
guidescan_guides <- guidescan_filtered %>% 
  group_by(id) %>% 
  mutate(guide_num = dplyr::row_number(), guide_id = paste(id, guide_num, sep = "_")) %>% 
  select(guide_id, id, gene, seq = gRNA.Seq, eff_score = Cutting.efficiency, chr, gRNA_start, 
         gRNA_end, guide_strand = Strand) 

#Reads in non-targeting guides from RE Screen performed in Gjaltema et al. 2022
NT_guides <- read.delim("./input_files/GjaltemaScreen_library_file.txt", 
                        col.names = c("guide_id", "seq", "gene")) %>% 
  filter(gene == "NT") %>% 
  mutate(id = "NT", chr = NA, gRNA_start = NA, gRNA_end = NA, guide_strand = NA, eff_score = NA) %>% 
  select(guide_id, id, gene, seq, eff_score, chr, gRNA_start, gRNA_end, guide_strand)


#Merges the guidescan and NT guides
total_guides <- full_join(guidescan_guides, NT_guides) %>% 
  write_delim("./output_files/TFi_total_guides.txt",
              delim = "\t", col_names = TRUE)

library_file <- total_guides %>% 
  select(guide_id, seq, id) %>% 
  write_delim("./output_files/TFi_library_file.txt",
              delim = "\t", col_names = FALSE)

#Removes isoforms without guides
df_final <- read.delim("./data/TFi_isoforms_raw.txt")

genes_final <- df_final %>% 
  filter(id %in% guidescan_guides$id)

write_delim(genes_final, "./output_files/TFi_isoforms_clean.txt",
            delim = "\t", col_names = TRUE)

#Prints the reduced transcription units included in the screen
df_reduced <- genes_final %>% 
  group_by(gene, chromosome, strand, cluster, id) %>%
  summarize(start = min(TSS), end = max(TSS))

write_delim(df_reduced, "./output_files/TFi_transcription_units.txt",
            delim = "\t", col_names = TRUE)

#Gives file with full sequences
genscript_seq <- total_guides %>% 
  ungroup() %>% 
  transmute(sequences = paste("ATCTTGTGGAAAGGACGAAACACCG", seq, "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGT", sep = ""))  %>% 
  write_delim("./output_files/TFi_genscript_seqs.txt", delim = "\t",
              col_names = FALSE)

#Gives Bed files with guide target regions
target_guide_bed <- guidescan_guides %>% 
  ungroup() %>% 
  transmute(chr, gRNA_start, gRNA_end, name = guide_id, score = 1000, strand = guide_strand) %>% 
  write_delim("./output_files/TFi_guides.bed", delim = "\t",
              col_names = FALSE)
