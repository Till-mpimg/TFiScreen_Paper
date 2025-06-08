#Runs the Diffbind analysis for differential peak calling 
library(tidyverse)
library(DiffBind)
library(GenomicRanges)
library(rtracklayer)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Read in peak file paths and store in a dataframe
Peaks <- list.files("./output_files/", pattern = ".*sgNT.bed", full.names = TRUE)
peak_df <- data.frame(Peaks) %>%
  mutate(Factor = str_extract(Peaks, "[A-Za-z0-9]*(?=\\_sgNT.bed)"))

#Read in BAM file paths and combine with peak file dataframe
bamReads <- list.files("./output_files/", pattern = "_dedup.bam$", full.names = TRUE)
cnt_prep <- data.frame(bamReads) %>%
  mutate(SampleID = str_extract(bamReads, "[A-Za-z0-9-]*(?=\\_dedup.bam)")) %>% 
  mutate(SampleID = str_remove(SampleID, "XX_D2_"))


diffbind_prep <- cnt_prep %>% 
  separate(SampleID, c("Factor", "Treatment", "Replicate"), sep = "-") %>% 
  inner_join(peak_df) %>%
  mutate(Tissue = "XX")

#Create dataframe containing all comparisons
comps <- c("Otx2", "Zic3", "Nfrkb", "Foxd3")
treatA <- c(rep("sgNT", 5))
treatB <- c("Otx2", "Zic3", "Nfrkb", "Foxd3")
Factor <- c("H3K4me3", "H3K27ac", "H3K4me1")

comp_treatments <- data.frame(comps, treatA, treatB) 
comp_df <- list(comps = comps, Factor = Factor) %>%
  expand.grid() %>%
  mutate(id = paste(comps, Factor, sep = "_")) %>%
  left_join(comp_treatments)

#Create a list of dataframes containing the necessary data for all comparisons
f_subset <- function(x) {
  subset(diffbind_prep, Factor == x[2] & Treatment %in% c(x[4], x[5]))
}

diffbind_list <- apply(comp_df, 1, f_subset)
names(diffbind_list) <- comp_df$id

diffbind_fun <- function(a)
{
  filtered_list <- diffbind_list[sapply(names(diffbind_list), grepl, pattern = paste0(a))]
  
  #Create DBA objects
  diffbind_sheet <- lapply(filtered_list, function(x) dba(sampleSheet = x))

  #Count reads to create a binding matrix
  diffbind_count <- lapply(diffbind_sheet, function(x) dba.count(x, bParallel=FALSE, summits = FALSE))

  #Create contrast between conditions
  diffbind_contrast <- lapply(diffbind_count, function(x) dba.contrast(x, categories=DBA_TREATMENT, 
                                                                       minMembers = 2, bNot = FALSE))
  
  #Perform differential binding analysis
  diffbind_analysis <- lapply(diffbind_contrast, function(x) dba.analyze(x, DBA_ALL_METHODS))
  
  #Export dba objects as .txt files
  diffbind_tracks <- lapply(diffbind_analysis, function(x) dba.report(x, method=DBA_ALL_METHODS, contrast = 1, th=0.05, bUsePval = T, DBA_DATA_FRAME))
  diffbind_df <- lapply(diffbind_tracks, function(x) as.data.frame(x))
  sapply(names(diffbind_df), function (x) write_delim(diffbind_df[[x]], paste0("./output_files/", x,"_diffbind.bed"), 
                                                      delim = "\t", col_names = FALSE))
}

sapply(Factor, diffbind_fun)
