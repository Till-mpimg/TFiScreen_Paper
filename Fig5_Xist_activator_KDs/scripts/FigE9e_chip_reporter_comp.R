#Script to compare public chip data, with Xist activator KD CUT&Tag and reporter screens
library(tidyverse)
library(egg)
library(gridExtra)
library(GenomicRanges)
library(rtracklayer)
library(Rsubread)
library(ggpubr)

order <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93","RE95", "RE96", "RE97", "RE127")
prox <- c("RE57L", "RE57M", "RE57R", "RE58")

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Load files
reporter_screen <- read.delim("./input_files/REPi_mean_results.txt") %>% 
  dplyr::select(Gene, RE = reporter, zfc, fdr) %>% 
  mutate(sig = ifelse(fdr <= 0.2 & zfc >= 0, "act",
                      "no"), 
         location = ifelse(RE %in% prox, "proximal", "distal")) %>% 
  filter(Gene %in% c("Otx2_p2", "Zic3_p1"))

re_locations <- read.delim("./input_files/re_total.bed", 
                           header = FALSE, col.names = c("chr", "start", "end", "RE", "score", "strand", 
                                                         "thickStart", "thickEnd", "color")) %>%
  mutate(RE = str_remove(RE, "_")) %>% 
  filter(RE %in% order) %>% 
  dplyr::select(chr, start, end, RE)

cnt_results <- read_delim("./output_files/CnT_Xist_REs_LFC.txt",
                          delim = "\t") %>% 
  filter(mark == "H3K27ac" & line %in% c("sgOtx2", "sgZic3")) %>% 
  mutate(RE = str_remove(RE, "_"), Gene = case_when(
    line == "sgOtx2" ~ "Otx2_p2",
    line == "sgZic3" ~ "Zic3_p1"
  ))

#Write RE file into genomic ranges and do overlaps to peak files
re_gr <-GRanges(
  seqnames = re_locations$chr,       
  ranges = IRanges(
    start = re_locations$start,      
    end = re_locations$end           
  ),
  RE = re_locations$RE              
)

#
feature_df <- re_locations %>% 
  transmute(GeneID = RE, Chr = chr, Start = start, End = end, Strand = ".")

bam_dir <- "/project/agsgpa/Till/chip_revision/20250314_Chiptest/revision_bam/"

setwd(bam_dir)
temp = list.files(pattern=paste0(".*.bam$"), full.names = TRUE)
temp_paired <- temp[grep('ZIC3', temp)]
temp_single <- temp[grep('OTX2', temp)]

feature_counts_single <- featureCounts(temp_single, annot.ext = feature_df, nthreads = 3, 
                                       allowMultiOverlap = TRUE)

feature_counts_paired <- featureCounts(temp_paired, annot.ext = feature_df, nthreads = 3, isPairedEnd = TRUE,
                                       allowMultiOverlap = TRUE)

tf_vec <- rep(c("Otx2_p2", "Zic3_p1"), 11)

counts <- data.frame(feature_counts_single$counts, feature_counts_paired$counts) %>% 
  cbind(length = feature_counts_single$annotation$Length) %>% 
  rownames_to_column("RE") %>% 
  pivot_longer(-c(RE, length), names_to = "sample", values_to = "counts") %>% 
  cbind(Gene = tf_vec) %>% 
  mutate(norm_counts = counts / length * 1000)


reporter_df <- counts %>%
  full_join(reporter_screen) %>% 
  group_by(sample) %>% 
  mutate(zscore = scale(norm_counts))


#Do a heatmap
heat_df <- reporter_df %>% 
  left_join(cnt_results) %>%
  mutate(cnt_res = -lfc) %>% 
  pivot_longer(c(zfc, zscore, cnt_res), names_to = "assay", values_to = "value")

lapply(c("Zic3_p1", "Otx2_p2"), function(x){
  pheat_df <- heat_df %>% 
    ungroup() %>% 
    filter(Gene == x ) %>% 
    dplyr::select(RE, assay, value) %>% 
    pivot_wider(names_from = RE, values_from = value) %>% 
    column_to_rownames("assay")
  
  pheat_mat <- pheat_df %>% 
    as.matrix()
  
  mat_col <- colorRampPalette(c("#5BBFBF", "#FFFFFF", "#F7A81F"))(9)
  
  rep_heat <- pheatmap(t(pheat_mat), 
                       border_color = "#C6C6C6",
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       fontsize = 6, 
                       color = mat_col, 
                       breaks = seq(-3, 3.6666667, 0.6666667),
                       cellwidth = 0.2 / 0.0353, 
                       cellheight = 0.2 / 0.0353,
                       legend = TRUE)
  
  ggsave(paste0("./output_files/FigE9e_", x, "_chip_screen_heatmap.pdf"), rep_heat, dpi = 300,
         useDingbats=FALSE)     
})

  
  