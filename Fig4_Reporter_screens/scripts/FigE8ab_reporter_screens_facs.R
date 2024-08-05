#Script analyzes the FlowFISH data from the reporter screen
library(flowCore)
library(openCyto)
library(tidyverse)
library(egg)
library(gridExtra)
library(scales)
library(EnvStats)
library(ggcyto)
library(viridis)
filter <- dplyr::filter

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6, margin = margin(0.05,0,0.05,0, "cm")),
                  strip.background = element_blank(), legend.title = element_blank()))


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Specifies location of FCS files
run1_dir <- paste0(wd, "./input_files/REPi_fcs/Run1/")
run2_dir <- paste0(wd, "./input_files/REPi_fcs/Run2/")
run3_dir <- paste0(wd, "./input_files/REPi_fcs/Run3/")
run4_dir <- paste0(wd, "./input_files/REPi_fcs/Run4/")


#Reads flowset for all runs individually and combines them together
read_facs <- function(x) {
flowset <- read.flowSet(path = x)

setwd(paste0(wd, "./data/"))

ggcyto::autoplot(flowset, x="FSC-A", y="SSC-A")

rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.5e5, 2e5), "SSC-A"= c(0.35e5, 1.5e5))
BoundaryFrame <- Subset (flowset[[1]], rectGate) 

autoplot(flowset, x="FSC-A", y="SSC-A") + geom_gate(rectGate)



sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"))
singlets_flowset <- Subset (flowset, rectGate %&% sing_g)

autoplot(flowset, x="FSC-A", y="FSC-H") + geom_gate(sing_g)


extr_frame <- fsApply(singlets_flowset, exprs, simplify = FALSE)
extr_df <- lapply(extr_frame, data.frame)
flow_df <- data.frame(do.call(rbind,extr_df)) %>%
  rownames_to_column("sample") %>% 
  mutate(GFP = GFP...510.20...A)  %>%
  select(sample, GFP)
}

run1 <- read_facs(run1_dir)
run2 <- read_facs(run2_dir)
run3 <- read_facs(run3_dir)
run4 <- read_facs(run4_dir)

flow_df <- rbind(run1, run2, run3, run4)
  
flow_df$sample <- str_replace(flow_df$sample, "_0.*", "")
flow_df_sep <- flow_df %>% 
  separate(sample, c("run", "rep", "type")) %>% 
  filter(GFP >= 0)

#Creates help dataframes used to plot the negative controls
help_run1 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run1") %>% 
  mutate(SP307 = GFP, SP379 = GFP, SP453 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:2), names_to = "type", values_to = "GFP")

help_run2 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run2") %>% 
  mutate(SP456 = GFP, SP530 = GFP, SP531 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:2), names_to = "type", values_to = "GFP")

help_run3 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run3") %>% 
  mutate(SP437 = GFP, SP451 = GFP, SP452 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:2), names_to = "type", values_to = "GFP")

help_run4 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run4") %>% 
  mutate(SP458 = GFP, SP459 = GFP, SP460 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:2), names_to = "type", values_to = "GFP")
 
flow_df_help <- rbind(help_run1, help_run2, help_run3, help_run4)


#Plots all histograms together
setwd(wd)

hist_plot <- flow_df_sep %>% 
  filter(type != "neg") %>% 
  ggplot() +
  facet_wrap(run~type, scales = "free", nrow = 2) +
  geom_density(data = flow_df_help, aes(x = GFP, ..scaled..), fill = "grey", color = NA) +
  geom_density(aes(x = GFP, ..scaled.., color = rep)) +
  scale_x_log10(breaks = c(10, 1000, 100000), limits = c(10, 100000), 
                labels = trans_format("log10", math_format(10^.x))) +
  scale_color_viridis(option = "F", discrete = TRUE, end = 0.8, begin = 0.2) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

fix <- set_panel_size(hist_plot, height = unit(1, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8a_REPi_facs.pdf", fix, dpi = 300,
       useDingbats=FALSE)




#Plots the mean MFI values of all samples together
mfi_df <- flow_df_sep %>% 
  group_by(run, rep, type) %>% 
  dplyr::summarize(mfi = geoMean(GFP)) %>%
  pivot_wider(names_from = type, values_from = mfi) %>% 
  mutate_at(vars(-run, -rep), list(~(. - neg))) %>% 
  select(-neg) %>% 
  pivot_longer(-c(1:2), names_to = "type", values_to = "mfi") %>% 
  na.omit()

mfi_plot <- mfi_df %>% 
  ggplot(aes(x = type, y = log10(mfi), color = rep)) +
  geom_point() +
  stat_summary(fun = "mean", geom = "crossbar", lwd = 0.25, width = 0.5, color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_viridis(option = "F", discrete = TRUE, end = 0.8, begin = 0.2) +
  scale_y_continuous(expand = c(0, 0.25))

fix <- set_panel_size(mfi_plot, height = unit(3, "cm"), width = unit(3.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8b_REPi_facs_mfi.pdf", fix, dpi = 300,
       useDingbats=FALSE)



#Returns 90th and 10th precentiles and MFI
neg_df <- flow_df_sep %>% 
  filter(type == "neg") %>% 
  group_by(run, rep, type) %>% 
  dplyr::summarize(mfi = geoMean(GFP)) %>%
  pivot_wider(names_from = type, values_from = mfi) 

perc_df <- flow_df_sep  %>% 
  filter(type != "neg") %>% 
  left_join(neg_df) %>% 
  mutate(GFP = GFP - neg) %>% 
  group_by(run, rep, type) %>% 
  summarize(low_gate = quantile(GFP, 0.1), high_gate = quantile(GFP, 0.9)) %>% 
  left_join(mfi_df) %>%  
  mutate(type = recode(type, `SP307` = "noRE", `SP379` = "RE57M", `SP437` = "RE97", `SP451` = "RE93",
                     `SP452` = "RE95", `SP453` = "RE96", `SP456` = "RE58", `SP458` = "RE61", `SP459` = "RE85",
                     `SP460` = "RE127", `SP530` = "RE57L", `SP531` = "RE57R"))

write_delim(perc_df, "./output_files/REPi_gates_mfi.txt", delim = "\t")




#Plots histograms for scheme

gates_df <- flow_df_sep  %>% 
  filter(type != "neg") %>% 
  group_by(run, rep, type) %>% 
  summarize(low_gate = quantile(GFP, 0.1), high_gate = quantile(GFP, 0.9))


hist_plot <- flow_df_sep %>% 
  filter(type != "neg") %>% 
  left_join(gates_df) %>% 
  filter(rep == "R1") %>% 
  ggplot() +
  facet_wrap(run~type, scales = "free") +
  geom_density(data = flow_df_help, aes(x = GFP, ..scaled..), fill = "grey", color = NA, alpha = 0.8) +
  geom_density(aes(x = GFP, ..scaled..), fill = "#72B854", color = "NA", alpha = 0.8) +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000), limits = c(100, 100000),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        axis.text.x = element_blank()) +
  geom_vline(aes(xintercept = low_gate)) +
  geom_vline(aes(xintercept = high_gate))

fix <- set_panel_size(hist_plot, height = unit(0.75, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4d_REPi_scheme.pdf", fix, dpi = 300,
       useDingbats=FALSE)
