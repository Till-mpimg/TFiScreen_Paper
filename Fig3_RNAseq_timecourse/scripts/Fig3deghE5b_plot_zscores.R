#Script to plot screen Tfs expression dynamics
library(tidyverse)
library(stats)
library(egg)
library(gridExtra)
library(EnvStats)

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


#Reads in data from TX timecourse and TFi Screen
clusters <- read.delim("./input_files/Timecourse_TF_clusters.txt")
cpm_tx <- read.delim("./input_files/CPM_RNA_timecourse.txt")
sig_df <- read.delim("./input_files/TFiLib_comp.txt")
tf_list <- read.delim("./input_files/Mus_musculus_TF.txt")$Symbol #From TF database
lib_file <- read.delim("./input_files/TFi_isoforms_clean.txt")
mle <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>%
  dplyr::rename(tu = "Gene") %>% 
  separate(tu, c("gene", "prom"), remove = FALSE)

TFi_TFs <- lib_file %>% 
  filter(gene %in% tf_list) %>%
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

mle_TFs <- mle %>% 
  filter(gene %in% TFi_TFs) %>% 
  group_by(gene) %>% 
  slice_max(abs(High.beta), n =1)


#Creates list of activators, repressors and others
act <- sig_df %>% 
  filter(sig != "ns" & TFi.beta < 0) %>%
  filter(gene %in% TFi_TFs) %>% 
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

rep <- sig_df %>% 
  filter(sig != "ns" & TFi.beta > 0) %>% 
  filter(gene %in% TFi_TFs) %>% 
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

rest <- TFi_TFs[!TFi_TFs %in% c(act, rep)] 

#Creates table with expression and only focuses on TFs and Xist
cpm_long <- cpm_tx %>%  
  filter(gene %in% c(act, rep, rest, "Xist")) %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("sex", "rep", "day"))

scale_cpm <- cpm_long %>% 
  filter(gene != "Xist") %>% 
  group_by(gene_id) %>% 
  mutate(zscore = as.numeric(scale(cpm)))

#Assign clusters back to genes
clust_df <- scale_cpm %>% 
  left_join(clusters) %>% 
  mutate(type = "ns")
clust_df[clust_df$gene %in% rep,]$type = "rep"
clust_df[clust_df$gene %in% act,]$type = "act"

#Export file
out_df <- clust_df %>% 
  select(-cpm) %>% 
  unite("sample", c(sex, rep, day)) %>% 
  pivot_wider(names_from = sample, values_from = zscore) %>% 
  mutate(cluster = fct_recode(cluster, transient_1 = "form2", transient_2 = "form1", transient_3 = "prime2", 
                              commited = "prime1", naive = "naive"))

write_delim(out_df, "./output_files/TFi_zscores_RNA_timecourse.txt", delim = "\t")

clust_df %>% select(gene, cluster, type) %>% unique() %>% group_by(cluster, type) %>% tally()

cluster_line_df <- clust_df %>%
  ungroup() %>% 
  select(-gene_id) %>% 
  mutate(day = str_remove(day, "h"))

cluster_smooth <- cluster_line_df %>%
  ggplot(aes(y = zscore, x = as.numeric(day))) +
  facet_wrap(sex~cluster, nrow = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.5) +
  stat_summary(aes(group = gene), color = "#666666", alpha = 0.1, geom = "line", fun = "mean") +
  stat_summary(data = cluster_line_df[cluster_line_df$type=="act",], aes(group = gene), color = "#F7A71C", alpha = 0.5,
               geom = "line", fun = "mean") +
  geom_smooth(aes(group = sex), method = "loess", se = FALSE, color = "black", 
              linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), limits = c(0, 96), name = "Timepoint [in hours]") +
  scale_y_continuous(name = "Z-score [expression]")

fix <- set_panel_size(cluster_smooth, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3d_Timecourse_cluster_expression.pdf", fix, dpi = 300,
       useDingbats=FALSE)


cluster_smooth_rep <- cluster_line_df %>%
  ggplot(aes(y = zscore, x = as.numeric(day))) +
  facet_wrap(sex~cluster, nrow = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.5) +
  stat_summary(aes(group = gene), color = "#666666", alpha = 0.1, geom = "line", fun = "mean") +
  stat_summary(data = cluster_line_df[cluster_line_df$type=="rep",], aes(group = gene), color = "#58BEBF", alpha = 0.5,
               geom = "line", fun = "mean") +
  geom_smooth(aes(group = sex), method = "loess", se = FALSE, color = "black", 
              linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), limits = c(0, 96), name = "Timepoint [in hours]") +
  scale_y_continuous(name = "Z-score [expression]")

fix <- set_panel_size(cluster_smooth_rep, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5b_Timecourse_cluster_expression_rep.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Plot TFi Screen data per TF cluster
clus_mle <- mle_TFs %>% 
  left_join(unique(clust_df[c("gene", "cluster", "type")])) %>% 
  arrange(High.beta) %>% 
  group_by(cluster) %>% 
  mutate(rank = order(High.beta, decreasing=FALSE))

#Calculate pvalue with fishers exact text
fisher_df <- expand.grid(cluster = unique(clus_mle$cluster),  type = c("act", "rep"))

fisher_df$pval <- mapply(function(a, b) {
  fun_df <- clus_mle
  fun_df[fun_df$cluster != a,]$cluster <- "rest"
  fun_df[fun_df$type != b,]$type <- "rest"
  
  p_val <- fisher.test(fun_df$cluster, fun_df$type, alternative = "greater")$p.value
  
  return(p_val)
}, fisher_df$cluster, fisher_df$type)

fisher_df$n_cluster <- sapply(fisher_df$cluster, function(a) {
  length(clus_mle[clus_mle$cluster==a,]$gene)
})

mle_plot <- clus_mle %>% 
  ggplot(aes(x = rank, y = High.beta)) +
  facet_wrap(~cluster, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("#F7A71C", "#58BEBF")) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-0.75, 0, 0.75), name = "Beta Score [TFi Lib]") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "none") +
  geom_segment(data = clus_mle[clus_mle$High.wald.p.value <= 0.05,], 
               aes(xend = rank, yend = 0, color = type)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_line(aes(group = cluster), color = "black") +
  geom_text(data = fisher_df[fisher_df$type=="act",], 
            aes(label = paste0("p=", round(pval, 3)), y = 0.75, x = 0.2*n_cluster), color = "#F7A71C", size = 6/2.8) +
  geom_text(data = fisher_df[fisher_df$type=="rep",], 
            aes(label = paste0("p=", round(pval, 3)), y = -0.75, x = 0.8*n_cluster), color = "#58BEBF", size = 6/2.8)


fix <- set_panel_size(mle_plot, height = unit(1, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3e_Timecourse_cluster_mle.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Check if clusters are biased for XX
sex_comp <- cluster_line_df %>% 
  filter(cluster %in% c("form1", "form2", "prime2")) %>%
  filter(type == "act")

cluster_sex <- sex_comp %>% 
  ggplot(aes(y = zscore, x = as.numeric(day))) +
  facet_wrap(~cluster, nrow = 1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.5) +
  geom_smooth(aes(group = sex, color = sex, fill = sex), method = "loess", se = TRUE, 
              linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), limits = c(0, 96), name = "Timepoint [in hours]") +
  scale_y_continuous(name = "Z-score [expression]", breaks = c(-1, 0, 1)) +
  scale_fill_manual(values = c("#CDCDCD", "#676767")) +
  scale_color_manual(values = c("#CDCDCD", "#676767"))

fix <- set_panel_size(cluster_sex, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3g_Timecourse_cluster_sex_onlyact.pdf", fix, dpi = 300,
       useDingbats=FALSE)



#Plot Screen categories
low_only <- c("Otx2", "Rbpj", "Foxd3", "Arnt", "Zfp217", "Arid1a", "Hmg20b", "Epas1", "Mbd3", "Zscan10", "Usf1")
high_low <- c("Oct4", "Zic3", "Nfrkb", "Nfe2l2", "Zfp518b", "Adnp", "Tox4", "Zfp207", "Sp1", "Pias1", "Myrf", 
              "Sall4", "Arid4a")

screen_comp <- cluster_line_df %>% 
  filter(gene %in% c(low_only, high_low)) %>%
  mutate(type = ifelse(gene %in% low_only, "low_only", "high_low"))

cluster_sex_screen <- screen_comp %>% 
  ggplot() +
  facet_wrap(~type, nrow = 1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.5) +
  geom_smooth(aes(group = sex, color = sex, fill = sex, y = zscore, x = as.numeric(day)), method = "loess", se = TRUE, 
              linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), limits = c(0, 96), name = "Timepoint [in hours]") +
  scale_y_continuous(name = "Z-score [expression]", breaks = c(-1, 0, 1)) +
  scale_fill_manual(values = c("#666666", "#F7A71C")) +
  scale_color_manual(values = c("#666666", "#F7A71C"))


fix <- set_panel_size(cluster_sex_screen, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3h_Timecourse_cluster_sex_screenRes.pdf", fix, dpi = 300,
       useDingbats=FALSE)