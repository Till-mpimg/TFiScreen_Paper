library(tidyverse)
library(egg)
library(gridExtra)

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
tpm_tx <- read.delim("./input_files/TPM_TX_timecourse_multifeat_pacini.txt")


sel_genes <- c("Pou5f1", "Zic3", "Nfrkb", "Otx2", "Foxd3")

tpm_long <- tpm_tx %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "tpm") %>% 
  separate(sample, c("replicate", "timepoint", "sex"))


#Plot selected factors
tpm_plot <- tpm_long %>% 
  mutate(timepoint = as.numeric(str_remove(timepoint, "h"))) %>% 
  filter(gene %in% sel_genes) %>% 
  ggplot(aes(x = timepoint, y = tpm, color = factor(sex, levels = c("XX", "XO")))) +
  facet_wrap(~gene, scales = "free") +
  stat_summary(fun = "mean", geom = "line") +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#F7A71C", "#666666"), name = "Sex") +
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), expand = c(0.1,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0.1, 0)) 


fix <- set_panel_size(tpm_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm")) 
grid.arrange(fix)
ggsave("./output_files/FigE5e_pacini_TPM.pdf", fix, dpi = 300,
       useDingbats=FALSE)
