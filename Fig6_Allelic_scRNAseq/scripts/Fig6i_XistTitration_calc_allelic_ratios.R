library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)
library(viridis)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]
setwd(wd)

xist_bin_col <- colorRampPalette(c("#ABABAB", "#93134D"))(6)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

cell_df <- read.delim("./output_files/XistTitration_filtered_cell_table.txt")

#Export suptable
sup_table <- cell_df %>% 
  dplyr::select(cell, sample, dTAG, rep, gene, count_xi, count_xa)
write_delim(sup_table, "./output_files/XistTitration_allelic_counts.txt", delim = "\t")

#Make totall allelic table
xist_info <- cell_df %>% 
  filter(gene == "Xist")

out_table <- cell_df %>%  
  filter(gene != "Xist") %>% 
  group_by(cell, sample, dTAG, rep, xi) %>% 
  summarize(count_xi = sum(count_xi), count_xa = sum(count_xa), n_genes = n()) %>% 
  mutate(allelic_counts = count_xi + count_xa, allelic_ratio = count_xi / (count_xi + count_xa))


#Plot Xist expression and allelic ratio per sample
sum_df <- cell_df %>%
  filter(gene != "Xist") %>% 
  group_by(Xist, cell, sample, dTAG) %>% 
  summarize(sum_xi = sum(count_xi), sum_xa = sum(count_xa)) %>% 
  mutate(chrX_ratio = sum_xi / (sum_xi + sum_xa))

test_df <- expand_grid(sample = unique(sum_df$sample), dTAG = unique(sum_df$dTAG)) %>% 
  filter(dTAG != "dTAG-500")

test_df$pval <- mapply(function(a,b) {
  wilcox.test(sum_df[sum_df$sample == a & sum_df$dTAG == "dTAG-500",]$chrX_ratio,
              sum_df[sum_df$sample == a & sum_df$dTAG == b,]$chrX_ratio)$p.value
}, test_df$sample, test_df$dTAG)

xci_plot <- sum_df %>% 
  ggplot(aes(x = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")), 
             y = chrX_ratio, 
             color = factor(dTAG, levels= c("dTAG-500", "dTAG-8", "dTAG-4", "dTAG-2", "dTAG-1", "dTAG-0")))) +
  facet_wrap(~sample) +
  geom_violin(fill = NA) +
  geom_jitter(size = 0.2, alpha = 0.2, width = 0.2) +
  stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.5, lwd = 0.25) +
  scale_color_manual(values = xist_bin_col) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_x_discrete(labels = c("500", "8", "4", "2", "1", "0"), name = "dTAG concentration [nM]")


allelic_fix <- set_panel_size(xci_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(allelic_fix)
ggsave("./output_files/Fig6i_XistTitration_chrX_ratio_dTAG.pdf", allelic_fix, 
       useDingbats=FALSE)