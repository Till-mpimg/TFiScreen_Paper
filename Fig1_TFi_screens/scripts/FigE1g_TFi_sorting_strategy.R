#This script plots the raw plots used for the TFi sorting strategy figure E1g

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
sample_flowset <- flowset[c("R2_TFi.fcs", "R2_2i.fcs")]

autoplot(sample_flowset, x="FSC-A", y="SSC-A")

#Set gate to remove dead cells and cell waste
rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.35e5, 1.5e5), "SSC-A"= c(0.25e5, 1.6e5))
BoundaryFrame <- Subset (sample_flowset[[1]], rectGate) 


autoplot(sample_flowset, x="FSC-A", y="SSC-A") + geom_gate(rectGate)

#Plots example for gating live cells
ssc_plot <- ggcyto(sample_flowset, aes(x = "SSC-A", y = "FSC-A"))+
  geom_hex(bins = 150) +
  ggcyto_par_set(limits = list(x = c(0, 250000), y= c(0, 250000))) +
  theme(axis.title.x = element_text(size = 6))+ 
  geom_gate(rectGate, colour = "darkred")


fix <- set_panel_size(as.ggplot(ssc_plot), height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1g_TFi_Sorting_a.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Set gate to remove doublets
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"))
singlets_flowset <- Subset (sample_flowset, rectGate %&% sing_g)

autoplot(sample_flowset, x="FSC-A", y="FSC-H") + geom_gate(sing_g)

#Plots example for singlet gating
sing_plot <- ggcyto(Subset (sample_flowset, rectGate), aes(x = "FSC-A", y = "FSC-H"))+
  geom_hex(bins = 150) +
  geom_gate(sing_g, colour = "black") +
  ggcyto_par_set(limits = list(x = c(0, 150000), y= c(0, 150000))) + 
  theme(axis.title.x = element_text(size = 6)) 


fix <- set_panel_size(as.ggplot(sing_plot), height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1g_TFi_Sorting_b.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Transform flowset back into dataframe
extr_frame <- fsApply(singlets_flowset, exprs, simplify = FALSE)

extr_df <- lapply(extr_frame, data.frame)
flow_df <- data.frame(do.call(rbind,extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub(".fcs.*", "", sample), Xist = APC...670.14...A) %>% 
  filter(Xist > 0) %>% 
  separate(sample, c("rep", "sample"), sep = "_")

quantile <- flow_df %>% 
  filter(sample == "2i") %>% 
  summarize(thresh =  quantile(Xist, 0.99)) %>% 
  unlist(use.names = FALSE)


#Plot screen data and controls by replicate as a histogram
hist_2i <- flow_df %>% 
  filter(sample == "2i") %>%  
  ggplot() +
  geom_density(aes(x = Xist, ..scaled..), fill = "#C6C6C6", color = NA) +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000), limits = c(50, 100000), name = "Xist RNA [Alexa 647, a.u.]",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.3), name = "Cell Density") +
  geom_vline(aes(xintercept = quantile), color = "darkred")

fix <- set_panel_size(hist_2i, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1g_TFi_Sorting_2i.pdf", fix, dpi = 300,
       useDingbats=FALSE)

neg_gate <- flow_df %>% 
  filter(sample == "TFi") %>% 
  summarize(thresh =  quantile(Xist, 0.15)) %>% 
  unlist(use.names = FALSE)

xist_vals <- flow_df %>% filter(sample == "TFi") %>% pull(Xist)
cdf <- ecdf(xist_vals)
thresh_quant <- cdf(quantile)
low_gate <- quantile(xist_vals, thresh_quant + 0.15)
high_gate <- quantile(xist_vals, 0.85)

#Plot TFi sorting
hist_TFi <- flow_df %>% 
  filter(sample == "TFi") %>%  
  ggplot() +
  geom_density(aes(x = Xist, ..scaled..), color = "black") +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000), limits = c(50, 100000), name = "Xist RNA [Alexa 647, a.u.]",
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.3), name = "Cell Density") +
  geom_vline(aes(xintercept = neg_gate), color = "darkred") +
  geom_vline(aes(xintercept = quantile), color = "darkred") +
  geom_vline(aes(xintercept = low_gate), color = "darkred") +
  geom_vline(aes(xintercept = high_gate), color = "darkred") 
  

fix <- set_panel_size(hist_TFi, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1g_TFi_Sorting_TFi.pdf", fix, dpi = 300,
       useDingbats=FALSE)
