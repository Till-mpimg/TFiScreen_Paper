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
xist_bin_col <- colorRampPalette(c("#ABABAB", "#93134D"))(6)


allelic_counts <- read.delim("./output_files/XistTitration_AllelicCounts_plusTotalXist.txt")
total_counts <- read.delim("./output_files/XistTitration_CountTable_plusMetadata_totalCounts.txt")


sup_table_total <- total_counts %>% 
  dplyr::select(cell = X, sample, Xist, replicate, total_counts, n_genes_by_counts, total_counts_mt) %>% 
  mutate(total_counts_mt = total_counts_mt / total_counts * 100)

write_delim(sup_table_total, "./output_files/XistTitration_totalcells_table.txt", delim = "\t")


xist_exp <- total_counts %>% 
  transmute(cell = X, sample, Xist, rep = replicate, xist_state = Xist_pos_neg) %>% 
  unique()

#Include gene names
gencode_xert <- "./input_files/GENCODE_vM25_plus_Xert_Linx.gtf"
gene_names <- read.gtf(gencode_xert) %>%
  filter(feature == "gene") %>% 
  mutate(position = ifelse(strand == "+", start, end)) %>% 
  dplyr::select(GeneID = gene_id, gene = gene_name, position) %>%
  unique()

#Load X chromosomal counts (chrX) 
cell_df <- allelic_counts %>% 
  dplyr::rename(GeneID = "gene") %>% 
  left_join(gene_names)

#Calculate %Xist+ cells per rep and sample
perc_df <- xist_exp %>% 
  separate(sample, c("sample", "dTAG"), sep = "_") %>% 
  group_by(sample, dTAG, rep, xist_state) %>% 
  tally() %>% 
  pivot_wider(names_from = xist_state, values_from = n) %>% 
  mutate(xist_perc = Xist_pos/(Xist_pos + Xist_neg))


test_df <- expand.grid(dTAG = c("dTAG-0", "dTAG-1", "dTAG-2", "dTAG-4", "dTAG-8"), 
                       sample = c("Xist-PROM", "Xist-RE"))

test_df$pval <- mapply(function(a, b) {
  t.test(perc_df[perc_df$sample==a & perc_df$dTAG==b,]$xist_perc, 
         perc_df[perc_df$sample==a & perc_df$dTAG=="dTAG-500",]$xist_perc, var.equal= TRUE)$p.value
}, test_df$sample, test_df$dTAG)

print(test_df)
#Plot Perc of xist_pos cells
xist_plot <- perc_df %>% 
  ggplot(aes(x = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")), 
             y = xist_perc, 
             color = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")))) +
  facet_wrap(~factor(sample, levels= c("Xist-PROM", "Xist-RE"))) +
  geom_point() +
  stat_summary(geom = "crossbar", fun = "mean", lwd = 0.25, width = 0.5, color = "black") +
  geom_text(data = test_df[test_df$pval <= 0.05,], aes(label = "*", y = 0.9), size = 8/2.8, color = "black") +
  scale_color_manual(values = xist_bin_col) +
  scale_y_continuous(limits = c(0, 1), name = "%Xist+ cells") +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("500", "8", "4", "2", "1", "0"), name = "dTAG concentration [nM]")



fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig6g_XistTitration_xistperc.pdf", fix, dpi = 300,
       useDingbats=FALSE)


# Gm14513 is filtered because all allelic reads map to ha1 (B6)
# Rps24-ps3 is filtered because all allelic reads map to ha2 (Cast)
fil_df <- cell_df %>% 
  left_join(xist_exp) %>% 
  filter(!gene %in% c("Gm14513", "Rps24-ps3"))


cell_n <- fil_df %>% dplyr::select(cell, sample, rep) %>% unique() %>% group_by(sample, rep) %>% tally()
cell_n
#Detect Xi via Xist expression
xist_df <- fil_df %>% 
  separate(sample, c("sample", "dTAG"), sep = "_") %>% 
  filter(gene == "Xist") %>% 
  mutate(xi = ifelse(count_ha1 > count_ha2, "ha1", "ha2")) %>% 
  mutate(total_allelic_xist = count_ha1 + count_ha2, 
         xist_ratio = ifelse(xi == "ha1", count_ha1 / (count_ha1 + count_ha2), count_ha2 / (count_ha1 + count_ha2))) %>% 
  mutate(xist_type = ifelse(xist_ratio <= 0.8, "biallelic", ifelse(xi == "ha1", "ha1", "ha2"))) %>% 
  dplyr::select(cell, sample, dTAG, Xist, xi, total_allelic_xist, xist_ratio, xist_type, count_ha1, count_ha2, rep)

join_xi <- xist_df %>% 
  dplyr::select(cell, rep, xi, total_allelic_xist, xist_type) %>% 
  unique()

#Calculate medians and foldchange
xist_df %>% group_by(sample, dTAG) %>% summarize(median = median(Xist)) %>% mutate(reverse_log = exp(median) - 1)


test_df$level_pval <- mapply(function(a, b) {
  wilcox.test(xist_df[xist_df$sample==a & xist_df$dTAG==b,]$Xist, 
              xist_df[xist_df$sample==a & xist_df$dTAG=="dTAG-500",]$Xist, var.equal= TRUE)$p.value
}, test_df$sample, test_df$dTAG)

#Xist plots
xist_violin <- xist_df %>% 
  ggplot(aes(x = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")), y = Xist/log(2), 
             color = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")))) +
  facet_wrap(~sample) +
  geom_violin(fill = NA) +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.5) +
  scale_color_manual(values = xist_bin_col) +
  stat_summary(geom = "crossbar", fun = "median", width = 0.5, color = "black") +
  scale_y_continuous(name = "Norm. counts [log2]", limits = c(0, NA)) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("500", "8", "4", "2", "1", "0"), name = "dTAG concentration [nM]") +
  geom_text(data = test_df[test_df$level_pval <= 0.05,], aes(label = "*", y = 7), size = 8/2.8, color = "black") 


fix_xist_violin <- set_panel_size(xist_violin, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix_xist_violin)
ggsave("./output_files/Fig6g_XistTitration_xistlvl.pdf", fix_xist_violin, 
       useDingbats=FALSE)

#Xist plots allelic
xist_scatter <- xist_df %>% 
  ggplot(aes(x = count_ha1, y = count_ha2, color = xist_type)) +
  facet_wrap(~sample) +
  geom_point(size = 0.5, alpha = 0.3) +
  scale_y_continuous(name = "Xist counts [Cast]", limits = c(0, 75)) +
  scale_x_continuous(name = "Xist counts [B6]", limits = c(0, 75)) +
  scale_color_manual(values = c("purple", "red", "blue", "#ABABAB"))


fix_scatter <- set_panel_size(xist_scatter, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix_scatter)
ggsave("./output_files/FigE10o_XistTitration_Xist_biallelic_scatter.pdf", fix_scatter, 
       useDingbats=FALSE)


#Filter cells with biallelic xist
filter_cells <- xist_df[xist_df$xist_type=="biallelic",]$cell


#Calculate allelic ratio
cor_df <- fil_df %>%
  separate(sample, c("sample", "dTAG"), sep = "_") %>% 
  inner_join(join_xi) %>% 
  na.omit() %>% 
  mutate(allelic_ratio = ifelse(xi == "ha1", count_ha1 / (count_ha1 + count_ha2), count_ha2 / (count_ha1 + count_ha2))) %>% 
  mutate(xi = as.factor(xi), xist_type = as.factor(xist_type))

#Export table
out_table <- cor_df %>%
  filter(gene != "Xist") %>% 
  group_by(cell, xi, sample, dTAG, total_allelic_xist, xist_type) %>% 
  summarize(sum_ha1 = sum(count_ha1), sum_ha2 = sum(count_ha2)) %>% 
  mutate(allelic_fraction = ifelse(xi == "ha1", sum_ha1 / (sum_ha1 + sum_ha2), sum_ha2 / (sum_ha1 + sum_ha2))) %>% 
  mutate(xi = fct_recode(xi, B6 = "ha1", Cast = "ha2"), xist_type = fct_recode(xist_type, B6 = "ha1", Cast = "ha2"))

write_delim(out_table, "/project/ag_schulz/Till/TFi_Paper/revision/output/scXist_Xistpos_table.txt", delim = "\t")

#Filters cells with allelic_ratio > 0.6 and biallelic Xist
fil_allelic <- cor_df %>%
  filter(gene != "Xist") %>% 
  group_by(cell, xi, sample, dTAG, total_allelic_xist) %>% 
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
  select(cell, sample, dTAG, rep, gene, count_xi, count_xa, allelic_ratio, position, Xist, xi) %>% 
  mutate(xi = fct_recode(xi, B6 = "ha1", Cast = "ha2"))

write_delim(fil_cor, "./output_files/XistTitration_filtered_cell_table.txt", delim = "\t")
