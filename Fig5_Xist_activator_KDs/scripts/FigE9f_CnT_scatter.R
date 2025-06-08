#Plots volcanoes vs the NTCs
library(tidyverse)
library(Rsubread)
library(GenomicRanges)
library(DESeq2)
library(ggrepel)
library(gridExtra)
library(egg)
library(Rgb)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

marks <- c("H3K4me3", "H3K4me1", "H3K27ac")
setwd("./output_files/")


#Count reads in NTC peaks and calculate CPM
count_reads <- function(x) {
  temp = list.files(pattern=paste0("XX_D2_", x, "_sg[A-Za-z0-9]*.bam$"), full.names = TRUE)
  peaks_bed <- read.delim(paste0(x, "_composite.bed"), header = FALSE)
  peaks_saf <- peaks_bed %>% 
    mutate(GeneID = paste0(V1, ":", V2, "-", V3), Strand = ".") %>% 
    select(GeneID, Chr = V1, Start = V2, End = V3, Strand)
  
  feature_counts <- featureCounts(temp, annot.ext = peaks_saf, isPairedEnd = TRUE,
                                  nthreads = 4, allowMultiOverlap = TRUE)
  counts <- data.frame(feature_counts$counts, feature_counts$annotation) %>%
    select(-Chr, -Start, -End, -Strand, -Length)
  names(counts) <- gsub(x = names(counts), pattern = "\\.", replacement = "_")
  names(counts) <- gsub(x = names(counts), pattern = "_bam", replacement = "")
  names(counts) <- gsub(x = names(counts), pattern = "XX_D2_", replacement = "")
  
  total_counts <- feature_counts$stat %>% 
    select(-Status) %>% 
    summarize_all(sum) %>% 
    pivot_longer(everything(), names_to = "sample", values_to = "total_counts") %>%
    mutate(sample = gsub("-", "_", str_extract(sample, '[A-Za-z0-9-]*')))
  
  cpm <- counts %>% 
    pivot_longer(-GeneID, names_to = "sample", values_to = 'counts') %>% 
    left_join(total_counts) %>% 
    mutate(cpm = (counts * 1000000)/ total_counts) %>% 
    dplyr::rename(peak_id = GeneID)
  
  
  return(cpm)
}

cpm_total <- do.call(rbind, lapply(marks, count_reads))

setwd(wd)

#add guides and contact genes
cpm_df <- cpm_total %>% 
  separate(sample, c("mark", "target")) %>% 
  select(-counts, -total_counts) %>% 
  pivot_wider(names_from = target, values_from = cpm) %>% 
  pivot_longer(c(sgOtx2, sgZic3, sgFoxd3, sgNfrkb), names_to = "target", values_to = "cpm")

  
#Read Diffbind files for pvalue
diffbind_files <- list.files(path = "./output_files", pattern = ".*diffbind.bed$", full.names = TRUE)

#Function to read all diffbind files
read_diffbind <- function(file) {
  bed_df <- read_tsv(file, col_names = FALSE) %>% 
    select(c(1,2,3,11))
  colnames(bed_df) <- c("chr", "start", "end", "FDR")
  bed_df$target <- str_extract(basename(file), "^[[:alnum:]]*")
  bed_df$mark <- str_extract(basename(file), "(?<=_)[[:alnum:]]*")
  
  out <- bed_df %>% 
    left_join(guide_df) %>% 
    mutate(start = start - 1) %>% 
    mutate(peak_id = paste0(chr, ":", start, "-", end)) %>% 
    select(-c(chr, start, end))
  return(out)
}

diffbind_df <- do.call(rbind, lapply(diffbind_files, read_diffbind)) %>% 
  na.omit()


combi_df <- full_join(cpm_df, diffbind_df) %>% 
  mutate(FDR = ifelse(is.na(FDR), 1, FDR)) %>% 
  filter(!peak_id == "chr15:99735630-99736116") %>% 
  filter(!str_detect(peak_id, "chrUn"))


combi_sum <- combi_df %>% 
  na.omit() %>% 
  group_by(mark, guide, target) %>% 
  summarize(n = n(), up = sum(cpm > sgNT & FDR <= 0.05), down = sum(cpm < sgNT & FDR <= 0.05))

#Plot scatter
H3K27ac_df <- combi_df %>% 
  filter(mark == "H3K27ac")


H3K4me3_df <- combi_df %>% 
  filter(mark == "H3K4me3")


H3K4me1_df <- combi_df %>% 
  filter(mark == "H3K4me1")

#Plot figures
H3K27ac_plot <- H3K27ac_df %>% 
  ggplot(aes(x = log10(LR15+0.01), y = log10(cpm+0.01))) +
  facet_wrap(~target) +
  geom_point(alpha = 0.1, size = 0.05, color = "#ABABAB") +
  geom_point(data = H3K27ac_df[H3K27ac_df$FDR <= 0.05,], color = "#d66b6b", size = 0.05, alpha = 0.2) +
  scale_x_continuous(limits = c(-1, 2)) +
  scale_y_continuous(limits = c(-1, 2)) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K27ac",], 
            aes(x = -0.3, y = 1.8, label = paste0("Up: ", up)), size = 5/2.8) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K27ac",], 
            aes(x = 1.3, y = -0.8, label = paste0("Down: ", down)), size = 5/2.8)

fix <- set_panel_size(H3K27ac_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9f_H3K27ac_scatter.pdf", dpi = 300, fix)



H3K4me3_plot <- H3K4me3_df %>% 
  ggplot(aes(x = log10(LR15+0.01), y = log10(cpm+0.01))) +
  facet_wrap(~target) +
  geom_point(alpha = 0.1, size = 0.05, color = "#ABABAB") +
  geom_point(data = H3K4me3_df[H3K4me3_df$FDR <= 0.05,], color = "#d66b6b", size = 0.05, alpha = 0.2) +
  scale_x_continuous(limits = c(-1, 2)) +
  scale_y_continuous(limits = c(-1, 2)) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K4me3",], 
            aes(x = -0.3, y = 1.8, label = paste0("Up: ", up)), size = 5/2.8) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K4me3",], 
            aes(x = 1.3, y = -0.8, label = paste0("Down: ", down)), size = 5/2.8)

fix <- set_panel_size(H3K4me3_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9f_H3K4me3_scatter.pdf", dpi = 300, fix)

H3K4me1_plot <- H3K4me1_df %>% 
  ggplot(aes(x = log10(LR15+0.01), y = log10(cpm+0.01))) +
  facet_wrap(~target) +
  geom_point(alpha = 0.1, size = 0.05, color = "#ABABAB") +
  geom_point(data = H3K4me1_df[H3K4me1_df$FDR <= 0.05,], color = "#d66b6b", size = 0.05, alpha = 0.2) +
  scale_x_continuous(limits = c(-1, 2)) +
  scale_y_continuous(limits = c(-1, 2)) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K4me1",], 
            aes(x = -0.3, y = 1.8, label = paste0("Up: ", up)), size = 5/2.8) +
  geom_text(data = combi_sum[combi_sum$mark=="H3K4me1",], 
            aes(x = 1.3, y = -0.8, label = paste0("Down: ", down)), size = 5/2.8)

fix <- set_panel_size(H3K4me1_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9f_H3K4me1_scatter.pdf", dpi = 300, fix)

