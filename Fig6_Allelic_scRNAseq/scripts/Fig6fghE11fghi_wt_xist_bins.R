library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)
library(viridis)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

xist_bin_col <- colorRampPalette(c("#ABABAB", "#93134D"))(4)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

cell_df <- read.delim("./output_files/filtered_cell_table.txt")

#Export suptable
sup_table <- cell_df %>% 
  select(cell, sample, rep, gene, count_xi, count_xa)
write_delim(sup_table, "./output_files/allelic_counts.txt", delim = "\t")

#Silencing statistics (Barros de Sousa 2019)
sil_statistics <- read.delim("./input_files/sousa_silencing.txt") %>% 
  select(gene = gene.name, half_time = half.time)

xo_exp <- read.delim("./input_files/CPM_RNA_timecourse.txt") %>% 
  select(gene, XO_R1_96h, XO_R2_96h, XO_R3_96h) %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "cpm") %>% 
  group_by(gene) %>% 
  summarize(xo_exp = mean(cpm))

#Bin wildtype cells according to Xist expression
bins <- cell_df %>% 
  filter(sample == "WT") %>% 
  select(cell, Xist) %>% 
  unique() %>% 
  mutate(xist_bin = factor(ntile(Xist, n = 4)))

bins %>% group_by(xist_bin) %>% tally()

bin_df <- cell_df %>% 
  filter(sample == "WT") %>% 
  left_join(bins)

#Plot Xist expression and allelic ratio per bin
sum_bin_df <- bin_df %>%
  filter(gene != "Xist") %>% 
  group_by(Xist, cell, xist_bin) %>% 
  summarize(sum_xi = sum(count_xi), sum_xa = sum(count_xa)) %>% 
  mutate(chrX_ratio = sum_xi / (sum_xi + sum_xa))


xci_plot <- sum_bin_df %>% 
  ggplot(aes(x = xist_bin, y = chrX_ratio, color = xist_bin)) +
  geom_violin(fill = NA) +
  geom_jitter(size = 0.2, alpha = 0.2, width = 0.2) +
  stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.5, lwd = 0.25) +
  scale_color_manual(values = xist_bin_col) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed")


allelic_fix <- set_panel_size(xci_plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(allelic_fix)
ggsave("./output_files/Fig6g_chrX_ratio_WT_bins.pdf", allelic_fix, 
       useDingbats=FALSE)


xist_plot <- sum_bin_df %>% 
  ggplot(aes(x = factor(xist_bin), y = Xist/log(2), color = factor(xist_bin))) +
  geom_violin(fill = NA) +
  geom_jitter(size = 0.2, alpha = 0.2, width = 0.2) +
  stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.5, lwd = 0.25) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = xist_bin_col) 

xist_fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(2.5, "cm"))
grid.arrange(xist_fix)
ggsave("./output_files/Fig6f_Xist_WT_bins.pdf", xist_fix, 
       useDingbats=FALSE)


pearson_wt <- cor(sum_bin_df$Xist, sum_bin_df$chrX_ratio)

scatter_bins <- sum_bin_df %>% 
  ggplot(aes(x = Xist, y = chrX_ratio, color = factor(xist_bin))) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", lwd = 0.25) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = xist_bin_col) +
  scale_y_continuous(limits = c(0, 0.65), breaks = c(0, 0.25, 0.5)) +
  scale_x_continuous(limits = c(0, NA)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", lwd = 0.25) +
  annotate(geom = "text", 
          x = 3.2, y = 0.6, label = paste0("r=", round(pearson_wt, 2)), color = "black", size = 6/2.8)
  
fix_scatter <- set_panel_size(scatter_bins, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_scatter)
ggsave("./output_files/WT_bins_scatter.pdf", fix_scatter, 
       useDingbats=FALSE)


#Summarize allelic ratio/Bin
sum_genes <- bin_df %>% 
  group_by(xist_bin, gene, position) %>% 
  summarize(n = n(), ratio = mean(allelic_ratio))

#Export Sup_tables
write_delim(sum_genes, "./output_files/allelic_genes.txt", delim = "\t")

#Filter genes with >=15 cells per category
filter_genes <- sum_genes %>%
  ungroup() %>%  
  group_by(gene) %>% 
  summarize(min = min(n)) %>% 
  filter(min >= 15) %>% 
  select(gene) %>% 
  unlist()

filter_sum <- sum_genes %>% 
  filter(gene %in% filter_genes)


#Pheatmap of individual genes
pheat_mat <- filter_sum  %>%
  select(-n, -position) %>%
  pivot_wider(names_from = xist_bin, values_from = ratio) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

dist_mat <- dist(pheat_mat)
hc <- hclust(dist_mat, method = "ward.D2")
plot(hc)
cluster <- cutree(hc, 5)
hc_df <- data.frame(hclust = as.factor(cluster)) %>% 
  mutate(hclust = fct_recode(hclust, weak_XCI = "1", early_XCI = "2", strong_XCI = "3", no_XCI = "4", Xist = "5"))

#Plot as heatmap
mat_col <- colorRampPalette(c("#ca0020", "#FFFFFF", "#0571b0"))(21)
ann_col <- substr(viridis(5, option = "viridis"), 1, 7)

target_heat <- pheatmap(pheat_mat[hc$order,], 
                        cluster_rows = FALSE, 
                        annotation_row = hc_df,
                        annotation_colors = list(cluster = c("1" = ann_col[1], "2" = ann_col[2], 
                                                             "3" = ann_col[3], "4" = ann_col[4],  "5" = ann_col[5])),
                        cluster_cols = FALSE,
                        border_color = NA,
                        fontsize = 6, 
                        color = mat_col, 
                        breaks = seq(0, 1, 0.04761905), 
                        cellwidth = 0.25 / 0.0353, 
                        cellheight = 0.2 / 0.0353)

ggsave("./output_files/Fig6h_WT_bins_hclust_heat.pdf", target_heat, 
       useDingbats=FALSE)

hc_merge <- hc_df %>% 
  rownames_to_column("gene")

cluster_sum <- filter_sum %>% 
  left_join(hc_merge)

#Export xist bin info to txt
write_delim(cluster_sum, "./output_files/filtered_genes_WT_bins.txt", delim = "\t")


#Create big dataframe with cluster and extra statistics
hc_order <- c("no_XCI", "weak_XCI", "strong_XCI", "early_XCI")

plot_df <-  cluster_sum %>% 
  left_join(sil_statistics) %>% 
  left_join(xo_exp) %>% 
  mutate(lin_dist = abs(position - 103483254)) %>% 
  filter(hclust != "Xist")


#Export allelic ratio info for Sup_tables
sup_df <- cluster_sum %>% 
  select(-position, -n) %>% 
  pivot_wider(names_from = xist_bin, values_from = ratio)
write_delim(sup_df, "./output_files/allelic_gene_cluster.txt", delim = "\t")




ratio_plot <- plot_df %>% 
  ggplot(aes(x = xist_bin, y = ratio, group = xist_bin, color = xist_bin)) +
  facet_wrap(~factor(hclust, levels = hc_order), nrow = 1) +
  geom_jitter(size = 0.1, color = "#676767", width = 0.3) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), name = "Allelic ratio [Xi/Xtotal]") +
  theme(legend.position = "none") +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_color_manual(values = xist_bin_col) 


fix_ratio_hc <- set_panel_size(ratio_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_ratio_hc)
ggsave("./output_files/Fig6h_WT_bins_hclust_ratios.pdf", fix_ratio_hc, 
       useDingbats=FALSE)

#Plot characteristics of the clusters
gene_df <- plot_df %>%
  ungroup() %>% 
  select(-xist_bin, -n, -ratio) %>% 
  unique()

sil_plot <- gene_df %>% 
  ggplot(aes(x = factor(hclust, levels = rev(hc_order)), 
             y = half_time, 
             group = factor(hclust, levels = rev(hc_order)))) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(size = 0.25, color = "#676767", width = 0.3) +
  scale_y_continuous(name = "Est. silencing halftime [days]", limits = c(0, NA)) +
  scale_x_discrete(name = "XCI cluster")

fix_sil <- set_panel_size(sil_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_sil)
ggsave("./output_files/FigE11i_HC_bins_sil.pdf", fix_sil, 
       useDingbats=FALSE)


dist_plot <- gene_df %>% 
  ggplot(aes(x = factor(hclust, levels = rev(hc_order)), 
             y = lin_dist / 1000000, 
             group = factor(hclust, levels = rev(hc_order)))) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(size = 0.25, color = "#676767", width = 0.3) +
  scale_y_continuous(name = "Distance to Xist TSS [Mb]", limits = c(0, NA)) +
  scale_x_discrete(name = "XCI cluster")

fix_dist <- set_panel_size(dist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_dist)
ggsave("./output_files/FigE11g_HC_bins_dist.pdf", fix_dist, 
       useDingbats=FALSE)



exp_plot <- gene_df %>% 
  ggplot(aes(x = factor(hclust, levels = rev(hc_order)), 
             y = log10(xo_exp + 1), 
             group = factor(hclust, levels = rev(hc_order)))) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  geom_jitter(size = 0.25, color = "#676767", width = 0.3) +
  scale_y_continuous(name = "XO Expression [CPM+1, log10]", limits = c(0, NA)) +
  scale_x_discrete(name = "XCI cluster")

fix_exp <- set_panel_size(exp_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_exp)
ggsave("./output_files/FigE11h_HC_bins_exp.pdf", fix_exp, 
       useDingbats=FALSE)
