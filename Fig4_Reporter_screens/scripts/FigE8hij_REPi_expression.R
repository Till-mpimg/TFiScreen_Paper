#Do analysis on the reporter screen output
library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)

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

zfc <- read.delim("./input_files/REPi_results.txt")
sig_TFi <-  read.delim("./input_files/TFiLib_comp.txt")
mle_TFi <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt")
cluster <- read.delim("./input_files/Timecourse_TF_clusters.txt")
expression <- read.delim("./input_files/CPM_RNA_timecourse.txt")

#Calculate Zscore of timecourse
zscore_cpm <- expression %>% 
  filter(gene %in% cluster$gene) %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("sex", "replicate", "timepoint")) %>%  
  mutate(timepoint = as.numeric(str_remove(timepoint, "h"))) %>% 
  group_by(gene) %>% 
  mutate(zscore = as.numeric(scale(cpm))) %>% 
  group_by(gene, sex, timepoint) %>% 
  summarize(zscore = mean(zscore))

#Filter for significant interactions
sig_interactions <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(zfc = mean(zfc)) %>% 
  select(gene = Gene, reporter, fdr, zfc) %>% 
  filter(fdr <= 0.2) %>% 
  unique() %>% 
  mutate(gene = str_remove(gene, "_p[1-3]"), type = ifelse(zfc >= 0, "act", "rep"), 
         cluster = ifelse(reporter %in% c("RE57L", "RE57M"), "proximal", 
                          ifelse(reporter %in% c("RE61", "RE85", "RE96"), "distal", "noRE-like"))) %>% 
  ungroup() %>% 
  select(gene, type, cluster) %>% 
  unique()

#Plot expression of repressors and activators per element
zscore_full <- left_join(zscore_cpm, sig_interactions) %>% 
  na.omit()

#Plot everything as heatmap
#get significant interactions by reporter
sig_reps <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(zfc = mean(zfc)) %>% 
  select(gene = Gene, reporter, fdr, zfc) %>% 
  filter(fdr <= 0.2) %>% 
  unique() %>% 
  mutate(gene = str_remove(gene, "_p[1-3]"), type = ifelse(zfc >= 0, "act", "rep"), 
         cluster = ifelse(reporter %in% c("RE57L", "RE57M"), "proximal", 
                          ifelse(reporter %in% c("RE61", "RE85", "RE96"), "distal", "noRE-like"))) %>% 
  ungroup() %>% 
  select(gene, type, cluster, reporter, zfc)

heat_df <- left_join(sig_reps, zscore_cpm)%>% 
  filter(type == "act") %>% 
  filter(sex == "XX")
  
#Order TF genes by time course
order_time <- heat_df %>% 
  select(gene, timepoint, zscore) %>% 
  unique() %>% 
  group_by(gene) %>% 
  slice_max(zscore, n = 1) %>% 
  arrange(timepoint, -zscore) %>% 
  select(gene) %>% 
  unlist(use.names = FALSE)

#Create heatmaps including gene per reporter screen
mat_col <- colorRampPalette(c("#FFFFFF", "#fcbba1", "#a50f15"))(10)

heat_fun <- function(a) {
  loop_mat <- heat_df %>% 
    filter(reporter == a) %>% 
    filter(timepoint <= 48) %>% 
    select(gene, timepoint, zscore) %>% 
    pivot_wider(names_from = timepoint, values_from = zscore) %>% 
    column_to_rownames("gene") %>% 
    as.matrix()
  
  dist <- dist(loop_mat)
  hclust <- hclust(dist, method = "ward.D")
  order <- order_time[order_time %in% row.names(loop_mat)]
  loop_order <- loop_mat[order,]
  
  heat <- pheatmap(loop_order, 
                       border_color = NA,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       fontsize = 6, 
                       color = mat_col, 
                       breaks = seq(-1, 2, 0.3),
                       cellwidth = 0.25 / 0.0353, 
                       cellheight = 0.05 / 0.0353,
                       treeheight_row = 0, 
                       show_rownames = FALSE,
                       legend = TRUE)
  
  
  ggsave(paste0("./output_files/E8h_", a, "_act_expression_heat.pdf"), heat, dpi = 300,
         useDingbats=FALSE)
  
}

lapply(c("RE57M", "RE57L", "RE61", "RE85", "RE96"), heat_fun)

#Plot as Boxplots
#Group by kmeans clusters (distal, proximal and noRE-like)
z_box <- zscore_full %>%
  filter(sex == "XX" & timepoint <= 48) %>% 
  ggplot(aes(x = factor(timepoint), y = zscore, 
             fill = factor(cluster, levels = c("proximal", "distal", "noRE-like")))) +
  facet_wrap(~type, nrow = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = c("#90134D", "#F7A71C", "#666666")) +
  scale_y_continuous(name = "Expression in XX [z-score]") +
  scale_x_discrete(name = "Timepoint after -2iL [hours]")

fix <- set_panel_size(z_box, height = unit(2, "cm"), width = unit(4, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8i_cluster_expression.pdf", fix, dpi = 300,
       useDingbats=FALSE)

test_df <- expand.grid(timepoint = c(0, 10, 16, 24, 30, 36, 48), cluster = c("proximal", "distal"), type = c("act", "rep"))
test_func <- function(a, b, c) {
  wilcox.test(zscore_full[zscore_full$timepoint==a & zscore_full$cluster==b & 
                            zscore_full$type==c & zscore_full$sex=="XX",]$zscore,
              zscore_full[zscore_full$timepoint==a & zscore_full$cluster==b & 
                            zscore_full$type==c & zscore_full$sex=="XO",]$zscore)$p.value
}

test_df$pval <- mapply(test_func, test_df$timepoint, test_df$cluster, test_df$type)

z_sex <- zscore_full %>%
  filter(timepoint <= 48 & cluster != "noRE-like") %>% 
  ggplot(aes(x = factor(timepoint), y = zscore)) +
  facet_wrap(cluster~type, nrow = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = factor(sex, levels = c("XX", "XO"))), 
               position = position_dodge(width = 0.8), outlier.size = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = c("#F7A71C", "#666666")) +
  scale_y_continuous(name = "Expression [z-score]") +
  scale_x_discrete(name = "Timepoint after -2iL [hours]") +
  geom_text(data = test_df[test_df$pval<= 0.05,], aes(y =2.2, label = "*"), size = 8/2.8, color = "black")

fix <- set_panel_size(z_sex, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8j_sex_expression.pdf", fix, dpi = 300,
       useDingbats=FALSE)
