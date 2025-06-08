#This script plots QC for TFiMini Screen
library(tidyverse)
library(egg)
library(gridExtra)
library(ggrepel)
library(viridis)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), 
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 6),  
                  strip.text = element_text(size = 6), strip.background = element_blank(),
                  legend.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Reads normalized count table
counts_complete <- read.delim("./input_files/TFiMini.count_normalized.txt")

counts <- counts_complete  %>%
  select(-c(1:2)) %>% 
  rename(PlasmidLib_PlasmidLib = PlasmidLib)


#Calculates distribution width of the populations excluding the NT guides
width_df <- counts_complete %>%
  filter(Gene != "NT") %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)

sample_levels = rev(c("PlasmidLib", "R1_Selected", "R1_2iSL", "R1_Negative",
                      "R1_High", "R2_Selected", "R2_2iSL", "R2_Negative", "R2_High"))

#Plots distribution widths
width_plot <- width_df %>%
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point() +
  labs(x = "Distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1, 3)

fix <- set_panel_size(width_plot, height = unit(2.5, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3d_TFiMini_log2_width.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Calculates distribution width of the fractions for NT guides only
ntc_width_df <- counts_complete %>%
  filter(Gene == "NT") %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)


#Plots NTC distribution widths
ntc_width_plot <- ntc_width_df %>%
  ggplot(aes(x = width, y = factor(sample, levels = sample_levels))) +
  geom_point() +
  labs(x = "NTC distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1, 3)


fix <- set_panel_size(ntc_width_plot, height = unit(2.5, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3d_TFiMini_NTC_log2_width.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Calculates Pearson Correlation for the normalized guide counts
counts_mat <- counts %>%
  as.matrix()

cor_mat <- round(cor(counts_mat),2)

cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sample2", values_to = "cor")

counts_long <- counts_complete %>%
  select(-PlasmidLib) %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("replicate", "fraction")) %>%
  pivot_wider(names_from = replicate, values_from = counts)

#dataframe to add the pearson correlations to the plot
cor_help <- cor_long %>% 
  separate(sample, c("repA", "fracA"), sep = "_") %>% 
  separate(sample2, c("repB", "fracB"), sep = "_") %>% 
  filter(fracA == fracB) %>% 
  filter(repA != repB) %>% 
  select(fraction = fracA, cor) %>% 
  unique()
  
#Plots the scatterplots for each screen population (as binned 2d density plot)
cor_plot <- counts_long %>%
  ggplot(aes(x = R1, y = R2)) +
  facet_wrap(~factor(fraction, levels = c("Selected", "2iSL", "Negative", "High")), nrow = 1) +
  geom_bin2d(bins = 60, aes(fill = after_stat(ndensity))) +
  geom_text(data = cor_help, aes(x= 3000, y = 12500, label = paste0("r=", cor)), size = 6/2.8) +
  labs(x = "Norm. counts (R1)", y = "Norm. counts (R2)") +
  scale_fill_viridis(option = "rocket") + 
  scale_y_continuous(breaks = c(0, 6000, 12000), limits = c(0, 15000)) +
  scale_x_continuous(breaks = c(0, 6000, 12000), limits = c(0, 15000))
  
fix <- set_panel_size(cor_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3e_TFiMini_correlation_density.pdf", fix, dpi = 300,
         useDingbats=FALSE)




