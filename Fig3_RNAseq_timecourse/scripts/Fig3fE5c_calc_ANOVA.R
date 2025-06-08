#Script to cluster TFs according to their expression profiles
library(tidyverse)
library(stats)
library(egg)
library(gridExtra)
library(EnvStats)
library(psych)

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
cpm_tx <- read.delim("./output_files/CPM_RNA_timecourse.txt")
sig_df <- read.delim("./input_files/TFiLib_comp.txt")
tf_list <- read.delim("./input_files/Mus_musculus_TF.txt")$Symbol #From TF database
lib_file <- read.delim("./input_files/TFi_isoforms_clean.txt")
clusters <- read.delim("./input_files/Timecourse_TF_clusters.txt")
mle <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>%
  dplyr::rename(tu = "Gene") %>% 
  separate(tu, c("gene", "prom"), remove = FALSE)

mle_low <- read.delim("./input_files/TFi_mle_LvN.gene_summary.txt") %>%
  dplyr::rename(tu = "Gene") %>% 
  separate(tu, c("gene", "prom"), remove = FALSE) %>% 
  group_by(gene) %>% 
  slice_max(abs(Low.beta), n =1) %>% 
  select(gene, Low_pval = Low.wald.p.value, Low_beta = Low.beta)


mle_high <- mle %>% 
  group_by(gene) %>% 
  slice_max(abs(High.beta), n =1) %>% 
  select(gene, High_pval = High.wald.p.value, High_beta = High.beta)

TFi_TFs <- lib_file %>% 
  filter(gene %in% tf_list) %>%
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)



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

test_df <- data.frame(gene = unique(cpm_long$gene))

test_df$pval_sex <- mapply(function(a) {
  test_aov <- cpm_long %>% 
    filter(gene == a)
  
  sum_aov <- summary(aov(cpm ~ sex * day, data = test_aov))
  p_value_sex <- sum_aov[[1]][["Pr(>F)"]][1]
  
  return(p_value_sex)
}, test_df$gene)

test_df$pval_diff <- mapply(function(a) {
  test_aov <- cpm_long %>% 
    filter(gene == a)
  
  sum_aov <- summary(aov(cpm ~ sex * day, data = test_aov))
  p_value_diff <- sum_aov[[1]][["Pr(>F)"]][2]
  
  return(p_value_diff)
}, test_df$gene)

test_df$fval_sex <- mapply(function(a) {
  test_aov <- cpm_long %>% 
    filter(gene == a)
  
  sum_aov <- summary(aov(cpm ~ sex * day, data = test_aov))
  p_value_sex <- sum_aov[[1]][["F value"]][1]
  
  return(p_value_sex)
}, test_df$gene)

test_df$fval_diff <- mapply(function(a) {
  test_aov <- cpm_long %>% 
    filter(gene == a)
  
  sum_aov <- summary(aov(cpm ~ sex * day, data = test_aov))
  p_value_diff <- sum_aov[[1]][["F value"]][2]
  
  return(p_value_diff)
}, test_df$gene)

test_df$fdr_sex <- p.adjust(test_df$pval_sex)
test_df$fdr_diff <- p.adjust(test_df$pval_diff)


means <- cpm_long %>% 
  group_by(gene, sex) %>% 
  summarize(cpm = mean(cpm)) %>% 
  pivot_wider(names_from = sex, values_from = cpm)

out_df <- test_df %>% 
  filter(gene != "Xist") %>% 
  left_join(means) %>% 
  mutate(lfc = log2(XX/XO), type = ifelse(gene %in% act, "act", ifelse(gene %in% rep, "rep", "rest"))) %>% 
  mutate(bias = ifelse(fdr_sex <= 0.05 & lfc >= 0, "XX", "no"), bias_diff = ifelse(fdr_diff <= 0.05, "bias", "unchanged")) 

write_delim(out_df, "./output_files/ANOVA_sex_tfs.txt", delim = "\t")

sum_sex <- out_df %>% 
  group_by(type, bias) %>% 
  tally() %>% 
  pivot_wider(names_from = bias, values_from = n) %>% 
  mutate(bias_perc = XX / (XX + no))


write_delim(sum_sex, "./output_files/ANOVA_sex_bias_summary.txt", delim = "\t")


sum_diff <- out_df %>% 
  group_by(type, bias_diff) %>% 
  tally() %>% 
  pivot_wider(names_from = bias_diff, values_from = n) %>% 
  mutate(bias_perc = bias / (bias + unchanged))


write_delim(sum_diff, "./output_files/ANOVA_diff_bias_summary.txt", delim = "\t")

#PLOT ANOVA results with screen results as dotplot
clus_levels <- c("naive", "form2", "form1", "prime2", "prime1")

dotplot_df <- out_df %>% 
  filter(type == "act") %>% 
  left_join(clusters) %>% 
  left_join(mle_high) %>% 
  left_join(mle_low) 

gene_order <- dotplot_df %>% 
  mutate(cluster = factor(cluster, levels = clus_levels)) %>% 
  arrange(cluster, High_beta) %>% 
  select(gene) %>% 
  unlist()
  

screen_df <- dotplot_df %>% 
  select(gene, cluster, High_pval, High_beta, Low_pval, Low_beta) %>% 
  pivot_longer(c(High_beta, Low_beta), names_to = "comp", values_to = "beta") %>% 
  pivot_longer(c(High_pval, Low_pval), names_to = "comp2", values_to = "pval") %>% 
  mutate(comp = str_extract(comp, "[A-Za-z]*"), comp2 = str_extract(comp2, "[A-Za-z]*")) %>% 
  filter(comp == comp2) %>% 
  select(-comp2) %>% 
  mutate(pval_bin = ifelse(pval > 0.2, ">0.2",
                           ifelse(pval > 0.05, "0.2", 
                                  ifelse(pval > 0.01, "0.05", "0.01"))))

mat_col <- colorRampPalette(c("#FFFFFF", "#F7A71C"))(10)
rna_col <- colorRampPalette(c("#FFFFFF", "#a50f15"))(10)

screen_dotplot <- screen_df %>% 
  ggplot(aes(x = comp, y = factor(gene, levels = rev(gene_order)), fill = -beta, size = pval_bin)) +
  geom_point(shape = 21) +
  scale_fill_stepsn(colors = mat_col, n.breaks = 10, limits = c(0, 1)) +
  scale_size_manual(values = c(">0.2" = 0.5, "0.2" = 1, "0.05" = 1.5, "0.01" = 2))
  

fix <- set_panel_size(screen_dotplot, height = unit(5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3f_Act_dotplot_screen.pdf", fix, dpi = 300,
       useDingbats=FALSE) 

anova_df <- dotplot_df %>% 
  select(gene, cluster, fval_sex, fval_diff, fdr_sex, fdr_diff) %>% 
  pivot_longer(c(fval_sex, fval_diff), names_to = "comp", values_to = "fval") %>% 
  pivot_longer(c(fdr_sex, fdr_diff), names_to = "comp2", values_to = "fdr") %>% 
  mutate(comp = str_extract(comp, "[A-Za-z]*$"), comp2 = str_extract(comp2, "[A-Za-z]*$")) %>% 
  filter(comp == comp2) %>% 
  select(-comp2) %>% 
  mutate(fdr_bin = ifelse(fdr > 0.05, ">0.05",
                           ifelse(fdr > 0.01, "0.05", 
                                  ifelse(fdr > 0.0001, "0.01", "0.0001"))))
  
  
anova_dotplot <- anova_df %>% 
  ggplot(aes(x = comp, y = factor(gene, levels = rev(gene_order)), fill = fval, size = fdr_bin)) +
  geom_point(shape = 21) +
  scale_fill_stepsn(colors = rna_col, n.breaks = 10, limits = c(0, 100)) +
  scale_size_manual(values = c(">0.05" = 0.5, "0.05" = 1, "0.01" = 1.5, "0.0001" = 2))


fix <- set_panel_size(anova_dotplot, height = unit(5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3f_Act_dotplot_anova.pdf", fix, dpi = 300,
       useDingbats=FALSE)   
  

#PLOT ANOVA results with screen results as dotplot for REPRESSORS
clus_levels <- c("naive", "form2", "form1", "prime2", "prime1")

dotplot_df <- out_df %>% 
  filter(type == "rep") %>% 
  left_join(clusters) %>% 
  left_join(mle_high) %>% 
  left_join(mle_low)

gene_order <- dotplot_df %>% 
  mutate(cluster = factor(cluster, levels = clus_levels)) %>% 
  arrange(cluster, -High_beta) %>% 
  select(gene) %>% 
  unlist()


screen_df <- dotplot_df %>% 
  select(gene, cluster, High_pval, High_beta, Low_pval, Low_beta) %>% 
  pivot_longer(c(High_beta, Low_beta), names_to = "comp", values_to = "beta") %>% 
  pivot_longer(c(High_pval, Low_pval), names_to = "comp2", values_to = "pval") %>% 
  mutate(comp = str_extract(comp, "[A-Za-z]*"), comp2 = str_extract(comp2, "[A-Za-z]*")) %>% 
  filter(comp == comp2) %>% 
  select(-comp2) %>% 
  mutate(pval_bin = ifelse(pval > 0.2, ">0.2",
                           ifelse(pval > 0.05, "0.2", 
                                  ifelse(pval > 0.01, "0.05", "0.01"))))

mat_col <- colorRampPalette(c("#FFFFFF", "#58BEBF"))(10)
rna_col <- colorRampPalette(c("#FFFFFF", "#a50f15"))(10)

screen_dotplot <- screen_df %>% 
  ggplot(aes(x = comp, y = factor(gene, levels = rev(gene_order)), fill = beta, size = pval_bin)) +
  geom_point(shape = 21) +
  scale_fill_stepsn(colors = mat_col, n.breaks = 10, limits = c(0, 1)) +
  scale_size_manual(values = c(">0.2" = 0.5, "0.2" = 1, "0.05" = 1.5, "0.01" = 2))


fix <- set_panel_size(screen_dotplot, height = unit(5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5c_dotplot_screen_REP.pdf", fix, dpi = 300,
       useDingbats=FALSE) 

anova_df <- dotplot_df %>% 
  select(gene, cluster, fval_sex, fval_diff, fdr_sex, fdr_diff) %>% 
  pivot_longer(c(fval_sex, fval_diff), names_to = "comp", values_to = "fval") %>% 
  pivot_longer(c(fdr_sex, fdr_diff), names_to = "comp2", values_to = "fdr") %>% 
  mutate(comp = str_extract(comp, "[A-Za-z]*$"), comp2 = str_extract(comp2, "[A-Za-z]*$")) %>% 
  filter(comp == comp2) %>% 
  select(-comp2) %>% 
  mutate(fdr_bin = ifelse(fdr > 0.05, ">0.05",
                          ifelse(fdr > 0.01, "0.05", 
                                 ifelse(fdr > 0.0001, "0.01", "0.0001"))))


anova_dotplot <- anova_df %>% 
  ggplot(aes(x = comp, y = factor(gene, levels = rev(gene_order)), fill = fval, size = fdr_bin)) +
  geom_point(shape = 21) +
  scale_fill_stepsn(colors = rna_col, n.breaks = 10, limits = c(0, 100)) +
  scale_size_manual(values = c(">0.05" = 0.5, "0.05" = 1, "0.01" = 1.5, "0.0001" = 2))


fix <- set_panel_size(anova_dotplot, height = unit(5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5c_dotplot_anova_REP.pdf", fix, dpi = 300,
       useDingbats=FALSE)   

