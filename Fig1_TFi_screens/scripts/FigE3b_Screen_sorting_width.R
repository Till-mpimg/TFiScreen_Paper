#Calculate Fold change during sorting and compare between replicates and screens
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

TFi_fcs <- "./input_files/TFi_FCS/"
TFiMini_fcs <- "./input_files/TFiMini_FCS/"

#Read FCS file as flowset
TFi_flowset <- read.flowSet(path = TFi_fcs, pattern="TFi")
TFiMini_flowset <- read.flowSet(path = TFiMini_fcs, pattern="TFi")

#Set gate to remove dead cells and cell waste
rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.35e5, 1.5e5), "SSC-A"= c(0.25e5, 1.6e5))
BoundaryFrame_TFi <- Subset (TFi_flowset[[1]], rectGate) 
BoundaryFrame_TFiMini <- Subset (TFiMini_flowset[[1]], rectGate) 

autoplot(TFiMini_flowset, x="FSC-A", y="SSC-A") + geom_gate(rectGate)

#Set gate to remove doublets
sing_g_TFi <- openCyto:::.singletGate(BoundaryFrame_TFi, channels = c("FSC-A", "FSC-H"))
singlets_flowset_TFi <- Subset (TFi_flowset, rectGate %&% sing_g)

sing_g_TFiMini <- openCyto:::.singletGate(BoundaryFrame_TFiMini, channels = c("FSC-A", "FSC-H"))
singlets_flowset_TFiMini <- Subset (TFiMini_flowset, rectGate %&% sing_g)

autoplot(flowset, x="FSC-A", y="FSC-H") + geom_gate(sing_g)

#Transform flowset back into dataframe
extr_frame_TFi <- fsApply(singlets_flowset_TFi, exprs, simplify = FALSE)
extr_frame_TFiMini <- fsApply(singlets_flowset_TFiMini, exprs, simplify = FALSE)

TFi_extr_df <- lapply(extr_frame_TFi, data.frame)
TFiMini_extr_df <- lapply(extr_frame_TFiMini, data.frame)

flow_df_TFi <- data.frame(do.call(rbind,TFi_extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub(".fcs.*", "", sample), Xist = APC...670.14...A) %>% 
  filter(Xist > 0) %>% 
  separate(sample, c("rep", "sample"), sep = "_")

flow_df_TFiMini <- data.frame(do.call(rbind,TFiMini_extr_df)) %>%
  rownames_to_column("sample") %>%
  transmute(sample = gsub(".fcs.*", "", sample), Xist = APC...670.14...A) %>% 
  filter(Xist > 0) %>% 
  separate(sample, c("rep", "sample"), sep = "_")

#Combine both screens and get quantiles
flow_df <- rbind(flow_df_TFi, flow_df_TFiMini) %>% 
  group_by(sample, rep) %>% 
  summarize(top_gate = quantile(Xist, 0.85), bottom_gate = quantile(Xist, 0.15)) %>% 
  mutate(lfc_gates = log10(top_gate) - log10(bottom_gate))

plot_gatewidth <- flow_df %>% 
  ggplot(aes(x = sample, color = sample, y = lfc_gates)) +
  geom_point() +
  scale_color_manual(values = c("orange", "brown")) +
  scale_y_continuous(limits = c(0, NA), expand = c(0.2, 0), name = "Xist detection width [Top - Bottom, log10]") +
  stat_summary(geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, color = "black")
  


fix <- set_panel_size(plot_gatewidth, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3b_Screen_Xist_gate_width.pdf", fix, dpi = 300,
       useDingbats=FALSE)

