#Do analysis on the reporter screen output
library(tidyverse)
library(egg)
library(gridExtra)

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
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")


#Creates list of activators, repressors and others
act <- sig_TFi %>% 
  filter(sig != "ns" & TFi.beta < 0) %>% 
  select(tu) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

rep <- sig_TFi %>% 
  filter(sig != "ns" & TFi.beta > 0) %>% 
  select(tu) %>% 
  unique() %>% 
  unlist(use.names = FALSE)


#Compares all interactions to mle
mean_zfc <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(mean_zfc = mean(zfc))

sig_zfc <- mean_zfc %>%  
  left_join(mle_TFi) %>% 
  mutate(type = ifelse(Gene %in% act, "act",
                       ifelse(Gene %in% rep, "rep", "ns")))

#Print numbers
n_sig <- sig_zfc %>% 
  filter(fdr <= 0.2) %>% 
  mutate(repi_type = ifelse(mean_zfc > 0, "act", "rep")) %>% 
  group_by(type, repi_type) %>% 
  tally()
n_sig

#Produce boxplot of all interactions
box_zfc <- sig_zfc %>%
  ggplot(aes(x = type, y = mean_zfc, color = type)) +
  geom_jitter(data = sig_zfc[sig_zfc$fdr >= 0.2,], size = 0.5, color = "#C6C6C6", alpha = 0.2, width = 0.3) +
  geom_jitter(data = sig_zfc[sig_zfc$fdr <= 0.2,], size = 0.5, alpha = 0.5, width = 0.3) +
  scale_y_continuous(limits = c(-3, 5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_color_manual(values = c("#F7A71C", "#666666", "#58BEBF"))
  
fix <- set_panel_size(box_zfc, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4e_zfc_boxplot.pdf", fix, dpi = 300,
       useDingbats=FALSE)
