#This script was used to produce plots describing TFi Library Design
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Reads targeted isoforms and output from Guidescan2 webtool
df_final <- read.delim("./output_files/TFi_isoforms_clean.txt")
guidescan_guides <- read.delim("./output_files/TFi_total_guides.txt") %>%
  filter(gene != "NT")

#Plots the number of TSSs per target promoter as a histogram
tss_per_TU <- df_final %>% 
  select(TSS, cluster, gene) %>% 
  unique() %>% 
  group_by(gene, cluster) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#C6C6C6", size = 0.1) +
  xlab("TSSs per TU")

fix <- set_panel_size(tss_per_TU, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1c_TSSs_per_TU.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Calculates the distance from the guide position to the targeted TSSs
guide_pos <- df_final %>%
  mutate(id = paste(gene, cluster, sep = "_")) %>% 
  ungroup() %>% 
  select(id, TSS, trans_strand = strand) %>% 
  inner_join(guidescan_guides, multiple = "all") %>% 
  mutate(dist_start = ifelse(trans_strand == "+", gRNA_start - TSS, TSS - gRNA_start)) %>% 
  mutate(dist_end = ifelse(trans_strand == "+", gRNA_end - TSS, TSS - gRNA_end)) %>% 
  mutate(dist = ifelse(abs(dist_start) < abs(dist_end), dist_start, dist_end )) %>% 
  select(id, TSS, guide_id, dist, trans_strand)/data/

#Plots a density plot of the guide distances to the TSSs
loc_plot <- guide_pos  %>% 
  ggplot() +
  geom_density(aes(x = dist, ..scaled..)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "#FAB336") +
  scale_y_continuous(name = "Guide Density [modal]", breaks = c(0, 0.5, 1)) +
  scale_x_continuous(name = "Distance to TSS [bp]", breaks = c(-500, 0, 500), limits = c(-700, 700))

fix <- set_panel_size(loc_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1d_guide_position.pdf", fix, dpi = 300,
       useDingbats=FALSE)
