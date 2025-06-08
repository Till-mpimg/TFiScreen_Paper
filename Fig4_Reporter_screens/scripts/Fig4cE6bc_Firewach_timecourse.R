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

#path to FCS files
run1a_dir <- paste0(wd, "./input_files/FCS_timecourse/Run1/load_a/")
run1b_dir <- paste0(wd, "./input_files/FCS_timecourse/Run1/load_b/")
run2_dir <- paste0(wd, "./input_files/FCS_timecourse/Run2/")
run3a_dir <- paste0(wd, "./input_files/FCS_timecourse/Run3/load_a/")
run3b_dir <- paste0(wd, "./input_files/FCS_timecourse/Run3/load_b/")

#Reads flowset for all runs individually and combines them together
read_facs <- function(x) {
  flowset <- read.flowSet(path = x)
  
  rectGate <-rectangleGate (filterId="Fluorescence Region", "FSC-A"= c(0.4e5, 1.8e5), "SSC-A"= c(0.2e5, 1.5e5))
  BoundaryFrame <- Subset (flowset[[1]], rectGate) 
  
  sqrcut <- matrix(c(35000,100000,250000,150000,50000,50000,150000,175000),ncol=2,nrow=4)
  colnames(sqrcut) <- c("FSC-A","FSC-H")
  
  sing_g <- polygonGate(BoundaryFrame, channels = c("FSC-A", "FSC-H"), boundaries = sqrcut)
  singlets_flowset <- Subset (flowset, rectGate %&% sing_g)
  
  extr_frame <- fsApply(singlets_flowset, exprs, simplify = FALSE)
  extr_df <- lapply(extr_frame, data.frame)
  flow_df <- data.frame(do.call(rbind,extr_df)) %>%
    rownames_to_column("sample") %>% 
    mutate(GFP = GFP...510.20...A)  %>%
    select(sample, GFP)
}

run1a <- read_facs(run1a_dir)
run1b <- read_facs(run1b_dir)
run2 <- read_facs(run2_dir)
run3a <- read_facs(run3a_dir)
run3b <- read_facs(run3b_dir)

#Combines flowsets into a single file
flow_df <- rbind(run1a, run1b, run2, run3a, run3b)

flow_df$sample <- str_replace(flow_df$sample, "_0.*", "")

#Renames samples to RE names
flow_df_sep <- flow_df %>% 
  separate(sample, c("run", "rep", "day", "type")) %>% 
  filter(GFP > 0) %>% 
  mutate(type = recode(type, SP307 = "noRE", SP379 = "RE57M", SP437 = "RE97", SP451 = "RE93", SP452 = "RE95", SP453 = "RE96",
         SP454 = "RE55", SP456 = "RE58", SP457 = "RE59", SP458 = "RE61", SP459 = "RE85", SP460 = "RE127",
         SP511 = "RE12", SP512 = "RE46", SP513 = "RE47", SP515 = "RE49", SP516 = "RE50",
         SP517 = "RE51", SP518 = "RE52", SP519 = "RE53", SP530 = "RE57L", SP531 = "RE57R"))

#Writes intensity values to output
write_delim(flow_df_sep, "./output_files/Firewach_Timecourse.txt", delim = "\t")

#Creates help dataframes used to plot the negative controls
help_run1 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run1") %>% 
  mutate(noRE = GFP, RE97 = GFP, RE93 = GFP, RE95 = GFP, RE96 = GFP, RE55 = GFP, RE58 = GFP, RE59 = GFP,
         RE61 = GFP, RE127 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:3), names_to = "type", values_to = "GFP")

help_run2 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run2") %>% 
  mutate(noRE = GFP, RE85 = GFP, RE12 = GFP, RE46 = GFP, RE47 = GFP, RE49 = GFP, 
         RE50 = GFP, RE51 = GFP, RE52 = GFP, RE53 = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:3), names_to = "type", values_to = "GFP")

help_run3 <- flow_df_sep %>%
  filter(type == "neg" & run == "Run3") %>% 
  mutate(noRE = GFP, RE57M = GFP, RE57L = GFP, RE57R = GFP) %>% 
  select(-GFP, -type) %>% 
  pivot_longer(-c(1:3), names_to = "type", values_to = "GFP")


flow_df_help <- rbind(help_run1, help_run2, help_run3) %>% 
  filter(day == "D0")


#Plots histograms of all replicates
#Only replicate 3 is shown in the paper
hist_plot <- flow_df_sep %>% 
  filter(type != "neg") %>%
  filter(!(type == "noRE" & run != "Run1")) %>% 
  ggplot() +
  facet_wrap(rep~type, nrow = 6) +
  geom_density(data = flow_df_help, 
               aes(x = GFP, ..scaled..), color = NA, fill = "grey") +
  geom_density(aes(x = GFP, ..scaled.., color = day), linewidth = 0.5) +
  scale_x_log10(breaks = c(10, 1000, 100000), limits = c(10, 100000),
                labels = trans_format("log10", math_format(10^.x)), name = "GFP [log10]") +
  scale_color_viridis(discrete = TRUE, end = 0.8) +
  scale_y_continuous(name = "Density [modal]", breaks = c(0, 0.5, 1.0)) +
  theme(legend.position = "top", legend.key.height = unit(0.125, "cm"), legend.key.width = unit(0.25, "cm"))

fix <- set_panel_size(hist_plot, height = unit(1, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE6c_Firewach_total_hist.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Plots direct comparison between empty and filled reporters
mfi_df <- flow_df_sep %>% 
  group_by(rep, day, type, run) %>% 
  dplyr::summarize(mfi = geoMean(GFP)) 

help_df <- mfi_df %>% 
  filter(type %in% c("noRE", "neg")) %>% 
  pivot_wider(names_from = type, values_from = mfi) %>% 
  mutate(noRE = noRE - neg)

#Writes relative MFI values to output
lfc_df <- mfi_df  %>% 
  filter(!type %in% c("noRE", "neg"))  %>%
  left_join(help_df) %>%  
  mutate(mfi = mfi - neg) %>% 
  mutate(lfc = log2(mfi/noRE)) %>% 
  select(-noRE, -neg)

write_delim(lfc_df, "./output_files/Firewach_Timecourse_relMFI.txt", delim = "\t")


mfi_plot <- lfc_df %>% 
  ggplot(aes(x = factor(day, levels = c("D0", "D1", "D2", "D3", "D4")), y = lfc)) +
  facet_wrap(~type, nrow = 2) +
  theme(legend.position = "none", panel.spacing = unit(0.05, "cm")) +
  geom_point(color = "grey", size = 0.5) +
  stat_summary(geom = "line", fun = "mean", aes(group = type), linewidth = 0.25) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0.5), name = "Rel. MFI to empty Reporter [log2]", breaks = c(-2, 0 , 2, 4)) +
  scale_x_discrete(name = "Day", expand = c(0, 0.5),
                   labels = c("0", "1", "2", "3", "4"))

fix <- set_panel_size(mfi_plot, height = unit(1, "cm"), width = unit(0.95, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3cE6b_Firewach_timecourse.pdf", fix, dpi = 300,
       useDingbats=FALSE)

