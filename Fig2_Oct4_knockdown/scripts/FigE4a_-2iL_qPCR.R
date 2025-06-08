#Script analyzing qPCRs from the validation of the TFi Screen
library(tidyverse)
library(scales)
library(egg)
library(EnvStats)
library(gridExtra)
library(ggh4x)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank(),
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


#Relative expression was calculated with ddCT method
rel_exp <- read.delim("./input_files/DIFF_relExp.txt")

qpcr <- rel_exp %>% 
  separate(Sample, c("rep", "day", "med"), sep = "_") %>% 
  pivot_longer(-c(rep, med, day), names_to = "gene", values_to = "rel_exp") %>% 
  na.omit()


#plot relative expression
y_list <- list(scale_y_continuous(expand = c(0,1), breaks = c(-15, -11, -7)),
               scale_y_continuous(expand = c(0,1), breaks = c(-4, -2, -0)),
               scale_y_continuous(expand = c(0,1), breaks = c(-1, 1, 3), limits = c(-1, 3)),
               scale_y_continuous(expand = c(0,1), breaks = c(-7, -5, -3)),
               scale_y_continuous(expand = c(0,1), breaks = c(-17, -13, -9)),
               scale_y_continuous(expand = c(0,1), breaks = c(-10, -5, 0)))

plot <- qpcr %>% 
  ggplot(aes(x = factor(day, levels = c("d0", "d1", "d2", "d3", "d4")), y = rel_exp, fill = med, group = med)) +
  facet_wrap(~gene, scales = "free") +
  geom_point(position = position_dodge(width = 1), shape = 21, stroke = 0.01) +
  stat_summary(fun = "mean", geom = "crossbar", color = "black", position = position_dodge(width = 1),
               lwd = 0.25) +
  facetted_pos_scales(y = y_list) +
  scale_fill_manual(values = c("#666666", "#F7A81F", "#58BEBF")) +
  scale_color_manual(values = c("white", "black"))

fix <- set_panel_size(plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE4a_-2iL_qPCR.pdf", fix, dpi = 300,
       useDingbats=FALSE)
