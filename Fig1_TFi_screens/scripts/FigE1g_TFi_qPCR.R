#Script to analyze QC qPCR for the TFi Screen
library(tidyverse)
library(scales)
library(egg)
library(gridExtra)
library(readxl)
library(EnvStats)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Takes the output .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC and removes all samples that did not include all target genes (e.g. negative controls).
raw <- read_xls("./input_files/TFi_screen_qPCR.xls", sheet = "Results", skip = 40) %>%
  select("Sample Name", "Target Name", "CT") %>%
  dplyr::rename(Sample = "Sample Name", Gene = "Target Name") %>%
  filter(Sample != "NTC") %>%
  mutate(CT = ifelse(CT == "Undetermined", 40, CT)) %>%  #Changes well with undetermined values to a CT of 40 (highest assessed cycle num)
  transform(CT = as.numeric(CT)) %>%
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarize(CT = geoMean(CT)) %>%
  pivot_wider(names_from = Gene, values_from = CT) %>%
  na.omit()


#Calculate relative expression compared to Arpo/Rrm2
analysis <- raw %>%
  mutate(Cont = (Arpo + Rrm2) /2, Arpo = NULL, Rrm2 = NULL) %>%
  select(Sample, Cont, everything()) 

for( i in 3:length(analysis) ) {
  analysis[i] <- analysis[i] - analysis[2]
  analysis[i] <- 2^ - analysis[i]
  analysis[i] <- log2(analysis[i])
}

output <- analysis %>%
  select(-Cont) 


plot_df <- output %>% 
  separate(Sample, c("rep", "sample"), sep = "_") %>% 
  pivot_longer(-c(1:2), names_to = "Gene", values_to = "Rel_Exp")

#Plots relative expression across gene controls for Xist, Prdm14 and Nanog
qpcr_plot <- plot_df %>% 
  ggplot(aes(x = sample, y = Rel_Exp, color = sample)) +
  facet_wrap(~Gene, scales = "free_y", nrow = 1) +
  geom_point() + 
  stat_summary(geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, color = "black") +
  scale_y_continuous(expand = c(0,1), name = "Relative expression (log2)") +
  scale_color_manual(values = c("#666666", "#90134D", "#C6C6C6", "#FAB336"))

fix <- set_panel_size(qpcr_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1g_TFi_qPCR.pdf", fix, dpi = 300,
       useDingbats=FALSE)
