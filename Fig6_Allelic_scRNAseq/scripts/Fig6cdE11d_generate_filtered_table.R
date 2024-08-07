library(tidyverse)
library(egg)
library(gridExtra)
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

xist_exp <- read.delim("./input_files/Xist_exp.txt")

#Include gene names
gencode_xert <- "./input_files/GENCODE_vM25_plus_Xert_Linx.gtf"
gene_names <- read.gtf(gencode_xert) %>%
  filter(feature == "gene") %>% 
  mutate(position = ifelse(strand == "+", start, end)) %>% 
  select(GeneID = gene_id, gene = gene_name, position) %>%
  unique()

#Load X chromosomal counts (chrX) 
cell_table <- "./input_files/allelic_chrX_counts.txt"
cell_df <- read.delim(cell_table) %>% 
  dplyr::rename(GeneID = "gene") %>% 
  left_join(gene_names)

#Calculate %Xist+ cells per rep and sample
perc_df <- xist_exp %>% 
  mutate(xist_state = ifelse(Xist>0, "Xist_pos", "Xist_neg")) %>% 
  group_by(sample, rep, xist_state) %>% 
  tally() %>% 
  pivot_wider(names_from = xist_state, values_from = n) %>% 
  mutate(xist_perc = Xist_pos/(Xist_pos + Xist_neg))


test_df <- data.frame(cont = c("WT"), 
                      sample = c("dFtx-Xert"))

test_df$pval <- mapply(function(a, b) {
  t.test(perc_df[perc_df$sample==a,]$xist_perc, perc_df[perc_df$sample==b,]$xist_perc, var.equal= TRUE)$p.value
}, test_df$cont, test_df$sample)

#pval = 0.14
#Plot Perc of xist_pos cells
xist_plot <- perc_df %>% 
  ggplot(aes(x = factor(sample, levels= c("WT", "dFtx-Xert")), 
             y = xist_perc, 
             color = factor(sample, levels= c("WT", "dFtx-Xert")))) +
  geom_point() +
  stat_summary(geom = "crossbar", fun = "mean", lwd = 0.25, width = 0.5, color = "black") +
  geom_text(data = test_df[test_df$pval <= 0.05,], aes(label = "*", y = 0.9), size = 8/2.8, color = "black") +
  scale_color_manual(values = c("#676767", "#E72B76")) +
  scale_y_continuous(limits = c(0, 1), name = "%Xist+ cells") +
  theme(axis.title.x = element_blank(), legend.position = "none")



fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig6c_scRNAseq_xistperc.pdf", fix, dpi = 300,
       useDingbats=FALSE)


# Gm14513 is filtered because all allelic reads map to ha1 (B6)
# Rps24-ps3 is filtered because all allelic reads map to ha2 (Cast)
fil_df <- cell_df %>% 
  left_join(xist_exp) %>% 
  filter(!gene %in% c("Gm14513", "Rps24-ps3"))

fil_df %>% select(cell, sample, rep) %>% unique() %>% group_by(sample, rep) %>% tally()

#Detect Xi via Xist expression
xist_df <- fil_df %>% 
  filter(gene == "Xist") %>% 
  mutate(xi = ifelse(count_ha1 > count_ha2, "ha1", "ha2")) %>% 
  mutate(total_allelic_xist = count_ha1 + count_ha2, 
         xist_ratio = ifelse(xi == "ha1", count_ha1 / (count_ha1 + count_ha2), count_ha2 / (count_ha1 + count_ha2))) %>% 
  mutate(xist_type = ifelse(xist_ratio <= 0.8, "biallelic", ifelse(xi == "ha1", "ha1", "ha2"))) %>% 
  select(cell, sample, Xist, xi, total_allelic_xist, xist_ratio, xist_type, count_ha1, count_ha2)

join_xi <- xist_df %>% 
  select(cell, xi, total_allelic_xist, xist_type)

#Calculate medians and foldchange
xist_df %>% group_by(sample) %>% summarize(median = median(Xist)) %>% mutate(reverse_log = exp(median) - 1)

#Xist plots
xist_violin <- xist_df %>% 
  ggplot(aes(x = sample, y = Xist/log(2), 
             color = sample)) +
  geom_violin(fill = NA) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.5) +
  scale_color_manual(values = c("#676767", "#E72E77")) +
  stat_summary(geom = "crossbar", fun = "median", width = 0.5, color = "black") +
  scale_y_continuous(name = "Norm. counts [log2]", limits = c(0, NA)) +
  theme(legend.position = "none", axis.title.x = element_blank())


fix_xist_violin <- set_panel_size(xist_violin, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix_xist_violin)
ggsave("./output_files/Fig6d_Xist_WT_DEL.pdf", fix_xist_violin, 
       useDingbats=FALSE)

#Xist plots allelic
xist_scatter <- xist_df %>% 
  ggplot(aes(x = count_ha1, y = count_ha2, color = xist_type)) +
  geom_point(size = 0.5, alpha = 0.3) +
  scale_y_continuous(name = "Xist counts [Cast]", limits = c(0, 32)) +
  scale_x_continuous(name = "Xist counts [B6]", limits = c(0, 32)) +
  scale_color_manual(values = c("purple", "red", "blue", "#ABABAB"))


fix_scatter <- set_panel_size(xist_scatter, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_scatter)
ggsave("./output_files/FigE11d_Xist_scatter.pdf", fix_scatter, 
       useDingbats=FALSE)


#Filter cells with biallelic xist
filter_cells <- xist_df[xist_df$xist_type=="biallelic",]$cell


#Calculate allelic ratio
cor_df <- fil_df %>%
  left_join(join_xi) %>% 
  na.omit() %>% 
  mutate(allelic_ratio = ifelse(xi == "ha1", count_ha1 / (count_ha1 + count_ha2), count_ha2 / (count_ha1 + count_ha2)))

#Export table
out_table <- cor_df %>%
  filter(gene != "Xist") %>% 
  group_by(cell, xi, sample, total_allelic_xist, xist_type) %>% 
  summarize(sum_ha1 = sum(count_ha1), sum_ha2 = sum(count_ha2)) %>% 
  mutate(allelic_fraction = ifelse(xi == "ha1", sum_ha1 / (sum_ha1 + sum_ha2), sum_ha2 / (sum_ha1 + sum_ha2))) %>% 
  mutate(xi = fct_recode(xi, B6 = "ha1", Cast = "ha2"), xist_type = fct_recode(xist_type, B6 = "ha1", Cast = "ha2"))

write_delim(out_table, "./output_files/Xistpos_table.txt", delim = "\t")

#Filters cells with allelic_ratio > 0.6 and biallelic Xist
fil_allelic <- cor_df %>%
  filter(gene != "Xist") %>% 
  group_by(cell, xi, sample, total_allelic_xist) %>% 
  summarize(sum_ha1 = sum(count_ha1), sum_ha2 = sum(count_ha2)) %>% 
  mutate(allelic_ratio = ifelse(xi == "ha1", sum_ha1 / (sum_ha1 + sum_ha2), sum_ha2 / (sum_ha1 + sum_ha2)),
         allelic_count = sum_ha1 + sum_ha2) %>% 
  filter(allelic_ratio >= 0.6) %>% 
  ungroup() %>% 
  select(cell) %>% 
  unlist()

fil_cor <- cor_df %>%
  filter(!cell %in% fil_allelic) %>% 
  filter(!cell %in% filter_cells) %>% 
  mutate(count_xi = ifelse(xi == "ha1", count_ha1, count_ha2), count_xa = ifelse(xi == "ha1", count_ha2, count_ha1)) %>% 
  select(cell, sample, rep, gene, count_xi, count_xa, allelic_ratio, position, Xist, xi) %>% 
  mutate(xi = fct_recode(xi, B6 = "ha1", Cast = "ha2"))

write_delim(fil_cor, "./output_files/filtered_cell_table.txt", delim = "\t")
