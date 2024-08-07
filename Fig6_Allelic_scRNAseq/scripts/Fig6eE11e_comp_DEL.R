library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)

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


cell_df <- read.delim("./output_files/filtered_cell_table.txt")
gene_clusters <- read.delim("./output_files/filtered_genes_WT_bins.txt") %>% 
  select(gene, hclust) %>% 
  unique()


#Plot allelic ratio per sample
sum_df <- cell_df %>%
  filter(gene != "Xist") %>% 
  group_by(Xist, cell, sample) %>% 
  summarize(sum_xi = sum(count_xi), sum_xa = sum(count_xa)) %>% 
  mutate(chrX_ratio = sum_xi / (sum_xi + sum_xa))

#Calculate medians
sum_df %>% group_by(sample) %>% summarize(median(chrX_ratio))

xci_plot <- sum_df %>% 
  ggplot(aes(x = sample, y = chrX_ratio, color = sample)) +
  geom_violin(fill = NA) +
  geom_jitter(size = 0.2, alpha = 0.2, width = 0.2) +
  stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.5, lwd = 0.25) +
  scale_color_manual(values = c("#676767", "#E72E77")) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed")


allelic_fix <- set_panel_size(xci_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(allelic_fix)
ggsave("./output_files/Fig6e_chrX_ratio_WT_DEL.pdf", allelic_fix, 
       useDingbats=FALSE)



#Summarize allelic ratio/Sample
#FIlter for genes with at least 16 cells in the dFtx-Xert sample
sum_genes <- cell_df %>% 
  group_by(sample, gene) %>% 
  summarize(n = n(), ratio = mean(allelic_ratio))

fil_genes <- sum_genes %>% 
  filter(sample == "dFtx-Xert") %>% 
  filter(n >= 16)

#Filter genes from binning analysis
filter_sum <- sum_genes %>%
  left_join(gene_clusters) %>% 
  na.omit() %>% 
  filter(gene %in% fil_genes$gene)


ratio_plot <- filter_sum %>% 
  filter(gene != "Xist") %>% 
  ggplot(aes(x = sample, y = ratio, group = sample, color = sample)) +
  geom_jitter(size = 0.1, color = "#676767", width = 0.3) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_text(data = filter_sum[filter_sum$gene == "Rlim",], aes(label = gene), color = "black", size = 4/2.8) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), name = "Allelic ratio [Xi/Xtotal]") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#676767", "#E72E77")) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed")


fix_ratio_hc <- set_panel_size(ratio_plot, height = unit(1.5, "cm"), width = unit(1, "cm"))
grid.arrange(fix_ratio_hc)
ggsave("./output_files/FigE11e_WT_DEL_gene_ratios.pdf", fix_ratio_hc, 
       useDingbats=FALSE)

