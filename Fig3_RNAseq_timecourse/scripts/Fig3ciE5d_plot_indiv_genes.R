#Script to plot screen Tfs expression dynamics
library(tidyverse)
library(stats)
library(egg)
library(gridExtra)
library(ggh4x)
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


#Read deseq2 file and cpm
cpm_tx <- read.delim("./input_files/CPM_RNA_timecourse.txt")
deseq_df <- read.delim("./input_files/DEseq2_total.txt")
sig_df <- read.delim("./input_files/TFiLib_comp.txt")
tf_list <- read.delim("./input_files/Mus_musculus_TF.txt")$Symbol #From TF database
lib_file <- read.delim("./input_files/TFi_isoforms_clean.txt")

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


#Creating plots for the expression of individual TFs
gene_vec <- c("Xist", "Pou5f1", "Zic3", "Nfrkb", "Otx2", "Foxd3")

plot_df <- cpm_long %>%
  dplyr::rename(timepoint = "day") %>% 
  left_join(deseq_df) %>% 
  mutate(timepoint = as.numeric(str_remove(timepoint, "h")))

y_list <- list(scale_y_continuous(breaks = c(0, 20, 40), limits = c(0, NA), expand = c(0.1, 0), name = "Expression [CPM]"),
               scale_y_continuous(breaks = c(0, 10, 20), limits = c(0, NA), expand = c(0.1, 0)),
               scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, NA), expand = c(0.1, 0)),
               scale_y_continuous(breaks = c(0, 300, 600), limits = c(0, NA), expand = c(0.1, 0)),
               scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, NA), expand = c(0.1, 0)),
               scale_y_continuous(breaks = c(0, 200, 400), limits = c(0, NA), expand = c(0.1, 0)))


plot <-  plot_df %>% 
  filter(gene %in% gene_vec) %>% 
  ggplot(aes(x = timepoint, y = cpm, color = sex)) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  geom_point(size = 0.2) +
  geom_smooth(se = FALSE, linewidth = 0.5) +
  facetted_pos_scales(y = y_list) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), name = "Timepoint [in hours]") +
  scale_color_manual(values = c("#CDCDCD", "#676767")) +
  geom_text(data = plot_df[plot_df$padj <= 0.05 & plot_df$gene %in% gene_vec,], aes(label = "*", y = 5),
            size = 6/2.8, color = "black")

fix <- set_panel_size(plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3ci_Timecourse_CPM_TFs.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Creating plots for the expression of other regulators outside top activators
plot_df2 <- plot_df %>%
  filter(gene %in% act) %>% 
  filter(!gene %in% gene_vec)




act_plot <- plot_df2 %>% 
  ggplot(aes(x = timepoint, y = cpm, color = sex)) +
  facet_wrap(~gene, scales = "free_y", nrow = 5) +
  geom_point(size = 0.2) +
  geom_smooth(se = FALSE, linewidth = 0.5) +
  scale_y_continuous(limits = c(0, NA), expand = c(0.1, 0)) +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), name = "Timepoint [in hours]") +
  scale_color_manual(values = c("#CDCDCD", "#676767")) +
  geom_text(data = plot_df2[plot_df2$padj <= 0.05,], aes(label = "*", y = 5),
            size = 6/2.8, color = "black")


fix <- set_panel_size(act_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5d_Timecourse_CPM_act.pdf", fix, dpi = 300,
       useDingbats=FALSE)