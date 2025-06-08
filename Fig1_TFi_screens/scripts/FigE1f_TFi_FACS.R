library(flowCore)
library(openCyto)
library(tidyverse)
library(egg)
library(gridExtra)
library(scales)
library(EnvStats)
library(ggcyto)
filter <- dplyr::filter

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  legend.key.height = unit(0.1, "cm"), legend.key.width = unit(0.25, "cm"),
                  strip.background = element_blank(), legend.title = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

fcs_dir <- "./input_files/TFi_FCS/"

#Read FCS file as flowset
flowset <- read.flowSet(path = fcs_dir)

autoplot(flowset, x="FSC-A", y="SSC-A")

#Set gate to remove dead cells and cell waste
rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.35e5, 1.5e5), "SSC-A"= c(0.25e5, 1.6e5))
BoundaryFrame <- Subset (flowset[[1]], rectGate) 

autoplot(flowset, x="FSC-A", y="SSC-A") + geom_gate(rectGate)

#Set gate to remove doublets
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"))
singlets_flowset <- Subset (flowset, rectGate %&% sing_g)

autoplot(flowset, x="FSC-A", y="FSC-H") + geom_gate(sing_g)

#Transform flowset back into dataframe
extr_frame <- fsApply(singlets_flowset, exprs, simplify = FALSE)

extr_df <- lapply(extr_frame, data.frame)
flow_df <- data.frame(do.call(rbind,extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub(".fcs.*", "", sample), Xist = APC...670.14...A) %>% 
  filter(Xist > 0) %>% 
  separate(sample, c("rep", "sample"), sep = "_")

#Set df for noProbe sample to use for plotting
help_df <- flow_df %>% 
  filter(sample == "noProbe") %>% 
  mutate(`2i` = Xist, TFi = Xist, sgNT = Xist, sgXist = Xist) %>% 
  select(-Xist, -sample) %>% 
  pivot_longer(-1, values_to = "Xist", names_to = "sample")

#Plot screen data and controls by replicate as a histogram
hist_plot <- flow_df %>% 
  filter(sample != "noProbe") %>%  
  ggplot() +
  facet_wrap(~rep, nrow = 1) +
  geom_density(data = help_df, 
               aes(x = Xist, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = Xist, ..scaled.., color = sample)) +
  scale_x_log10(breaks = c(10, 1000, 100000), limits = c(10, 100000), name = "Xist RNA [Alexa 647, a.u.]",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), name = "Cell Density") +
  scale_color_manual(values = c("#666666", "#90134D", "#C6C6C6", "#FAB336"))

fix <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1f_TFi_FACS.pdf", fix, dpi = 300,
       useDingbats=FALSE)
