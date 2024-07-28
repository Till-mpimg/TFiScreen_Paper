#Script to plot FISH data (Control for TFi Screen)
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Read FISH data
fish <- read.delim("./input_files/Fig_E1/TFi_FISH.txt")


#Plot the data as a barplot
xist_plot <- fish %>% 
  select(sample, rep, xist_pos, xist_2) %>% 
  pivot_longer(c(3:4), names_to = "num", values_to = "perc") %>% 
  ggplot(aes(x = sample, y = perc, fill = factor(num, levels = c("xist_pos", "xist_2")))) +
  facet_wrap(~rep, nrow = 1) +
  geom_bar(stat = "identity", color = "black", width = 0.7, size = 0.1) +
  scale_y_continuous(limits = c(0, 100), name = "%Xist-positive cells") +
  scale_fill_manual(values = c("#C6C6C6", "#666666"))

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1f_TFi_FISH.pdf", fix, dpi = 300,
       useDingbats=FALSE)

