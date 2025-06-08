#Designs the library for the TFiMini and the reporter screens
library(tidyverse)
library(readxl)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


#Reads mle from High/Negative analysis from TFi screen
mle_HvN <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>% 
  mutate(rank = rank(High.beta), ties.method = "first")


#Reads guide counts from the TFi screen
counts_complete <- read.delim("./input_files/TFi_corr.count_normalized.txt") %>% 
  select(sgRNA, Gene, R1_Negative, R2_Negative, R1_High, R2_High) %>%
  mutate(mean_counts = (R1_Negative + R2_Negative + R1_High + R2_High) / 4) %>% 
  filter(mean_counts > 10) %>% 
  transmute(sgRNA, Gene, lfc = log2((R1_High + R2_High) / (R1_Negative + R2_Negative)))

#Reads in screen counts from RE screen (Gjaltema 2022) and calculates fold change High / Negative
counts_RE <- read.delim("./input_files/Gjaltema_REscreen_norm_counts.txt") %>% 
  select(sgRNA = ID, RE, R1_Negative, R2_Negative, R1_High, R2_High) %>%
  mutate(mean_counts = (R1_Negative + R2_Negative + R1_High + R2_High) / 4) %>% 
  filter(mean_counts > 10) %>% 
  transmute(sgRNA, RE, lfc = log2((R1_High + R2_High) / (R1_Negative + R2_Negative)))

#Reads TF list
tfs <- read.delim("./input_files/Mus_musculus_TF_AnimalTFdb.txt")

#Add published TF regulators back to the library, even if they did not score
pub_tf <- c("Yy1", "Ctcf", "Nanog", "Zfp42", "Pou5f1", "Sox2", "Myc", "Klf4", "Prdm14", "Esrrb")

#Selects the highest scoring promoter of each Gene that had at least a pvalue <= 0.2 in the TFi screen
#Removes targets from the xist locus, as they will be added back later from the RE screen data
sig_p <- mle_HvN %>% 
  separate(Gene, c("Gene", "prom"), sep = "_") %>% 
  filter(High.wald.p.value <= 0.2 | Gene %in% pub_tf) %>%
  filter(!Gene %in% c("Xist", "Tsix", "Jpx", "Ftx", "Xert", "Rlim")) %>% 
  group_by(Gene) %>% 
  slice_min(High.wald.p.value) %>% 
  mutate(reg = ifelse(High.beta < 0, "up", "down")) %>%
  unite("Gene", c(Gene, prom)) %>% 
  select(Gene, reg)


#Selects the top 8 guides in terms of absolute fold change from the TFi screen
sel_guides_tfs <- counts_complete %>% 
  filter(Gene != "NT") %>% 
  filter(Gene %in% sig_p$Gene) %>% 
  group_by(Gene) %>% 
  slice_max(abs(lfc), n = 8)

#Selects the 100 NTs with the lowest absolute fold change from the TFi screen
sel_guides_nt <- counts_complete %>% 
  filter(Gene == "NT") %>% 
  slice_min(abs(lfc), n = 100)

#Selects the top 10 guides for selected REs
RE_list <- c("RE12", "RE46", "RE47", "RE49", "RE50", "RE51", "RE53", "RE57", "RE58", "RE61", "RE85", 
             "RE93", "RE95", "RE96", "RE97", "RE127")

sel_guides_RE <- counts_RE %>% 
  filter(RE %in% RE_list) %>% 
  group_by(RE) %>% 
  slice_max(abs(lfc), n = 10)


#Returns all guide sequences and completes the library
TFi_guide_seq <- read.delim("./input_files/TFi_library_file.txt",
                            header = FALSE, col.names = c("sgRNA", "seq", "Gene"))

Elib_guides <- left_join(sel_guides_RE, Elib_guide_seq) %>% 
  select(sgRNA, Gene = RE, seq)

TFi_guides <- rbind(data.frame(sel_guides_tfs),data.frame(sel_guides_nt)) %>% 
  left_join(TFi_guide_seq) %>% 
  select(sgRNA, Gene, seq)


#Adds guides that target Fgf4 minimal promoter as positive control (designed with CHOPCHOP)
fgf4_guides <- read.csv("./input_files/fgf4_guides.csv") %>% 
  rownames_to_column("num") %>% 
  transmute(sgRNA = paste0("FIREWACh_", num), Gene = "FIREWACh", seq = gRNA.Seq)


#Completes library and returns library file + genscript list
total_guides <- rbind(data.frame(Elib_guides), TFi_guides, fgf4_guides)

library_file <- total_guides %>% 
  write_delim("./output_files/TFiMini_library_file.txt",
              delim = "\t", col_names = FALSE)


#Add overhangs for genscript order
genscript_seq <- total_guides %>% 
  ungroup() %>% 
  transmute(sequences = paste("ATCTTGTGGAAAGGACGAAACACCG", seq, "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTAATGGACATCTTATTCACAG", sep = ""))


write_delim("./output_files/TFiMini_genescript_seqs.txt", delim = "\t",
              col_names = FALSE)

