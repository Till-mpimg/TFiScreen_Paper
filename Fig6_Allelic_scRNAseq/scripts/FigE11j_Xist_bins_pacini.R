#Compare with clusters from Schwaemmle 2024 scrnaseq
library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)
library(Rgb)
library(viridis)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

gene_clusters <- read.delim("./input_files/filtered_genes_WT_bins.txt") %>% 
  select(gene, hclust) %>% 
  unique()

#Filter Gm14513 and Rps4x because all allelic reads map to the same allele (misassigned SNP?)
cell_df <- read.delim("./output_files/AllelicCounts_pacini_day3.txt") %>% 
  filter(!gene %in% c("Gm14513", "Rps4x"))

xist_df <- cell_df %>% 
  filter(gene == "Xist") %>% 
  filter(count_B6 != count_Cast) %>% 
  mutate(xist_allelic_counts = count_B6 + count_Cast, Xi = ifelse(count_B6 > count_Cast, "B6", "Cast")) %>% 
  filter(xist_allelic_counts > 5) %>% #Filter Xist allelic reads
  mutate(xist_ratio = ifelse(Xi == "B6", count_B6 / (xist_allelic_counts), count_Cast / (xist_allelic_counts))) %>% 
  select(cell, xist_allelic_counts, Xi, xist_ratio) %>%
  mutate(xist_type = ifelse(xist_ratio <= 0.8, "biallelic", ifelse(Xi == "B6", "B6", "Cast")))

#Combine Xist information into cell_df
comb_df <- cell_df %>% 
  filter(cell %in% xist_df$cell) %>% 
  left_join(xist_df)


#Calculate allelic ratio
cor_df <- comb_df %>%
  filter(xist_type != "biallelic") %>% 
  na.omit() %>% 
  mutate(allelic_ratio = ifelse(Xi == "B6", count_B6 / (count_B6 + count_Cast), count_Cast / (count_B6 + count_Cast))) %>% 
  mutate(Xist_bin = factor(ntile(Xist_norm, 4)))

sup_cells <- cor_df %>% 
  select(cell, Xist_norm, xist_allelic_counts, Xi, Xist_bin) %>% 
  unique() %>% 
  mutate(Xist_norm = Xist_norm / log(2))

write_delim(sup_cells, "./output_files/pacini_cell_table.txt", delim = "\t")

sup_counts <- cor_df %>% 
  select(cell, gene, count_B6, count_Cast, xist_type) %>% 
  filter(count_B6 + count_Cast >= 1)

write_delim(sup_counts, "./output_files/pacini_count_table.txt", delim = "\t")
  

#Summarize allelic ratio/gene
sum_genes <- cor_df %>% 
  na.omit() %>% 
  group_by(Xist_bin, gene) %>% 
  summarize(n = n(), ratio = mean(allelic_ratio))

#Export gene data
write_delim(sum_genes, "./output_files/pacini_gene_table.txt", delim = "\t")

#Filter genes in new scRNAseq
filter_sum <- sum_genes %>% 
  left_join(gene_clusters) %>% 
  na.omit() %>% 
  filter(gene != "Xist") %>% 
  filter(n >= 8)


ratio_plot <- filter_sum %>% 
    ggplot(aes(x = Xist_bin, y = ratio, group = Xist_bin, color = Xist_bin)) +
    facet_wrap(~factor(hclust, levels = hc_order), nrow = 1) +
    geom_jitter(size = 0.1, color = "#676767", width = 0.3) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), name = "PACINI Allelic ratio [Xi/Xtotal]") +
    theme(legend.position = "none") +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
    scale_color_manual(values = xist_bin_col) 
  
  
  fix_ratio_hc <- set_panel_size(ratio_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
  grid.arrange(fix_ratio_hc)
  ggsave("./output_files/FigE11j_PACINI_hclust_ratios.pdf", fix_ratio_hc, 
         useDingbats=FALSE)
  