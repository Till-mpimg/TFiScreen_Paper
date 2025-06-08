library(tidyverse)
library(egg)
library(gridExtra)
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





