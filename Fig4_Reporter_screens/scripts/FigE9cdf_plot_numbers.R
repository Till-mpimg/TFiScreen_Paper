#Do analysis on the reporter screen output
library(tidyverse)
library(egg)
library(gridExtra)
library(EnvStats)

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
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")
facs_results <- read.delim("./input_files/REPi_gates_mfi.txt")


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


#Plot per element
n_df <- mean_zfc %>% 
  filter(fdr <= 0.2) %>% 
  mutate(type = ifelse(mean_zfc <= 0, "rep", "act")) %>% 
  group_by(type, reporter) %>% 
  tally()

n_plot <- n_df %>% 
  ggplot(aes(x = factor(reporter, levels = level), y = n, fill = type)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values = c("#F7A71C", "#58BEBF")) +
  scale_y_continuous(name = "Number of interactions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank())

fix <- set_panel_size(n_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9c_n_interactions.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Plots the number of detected interactions vs the distance between the high/low gates
mean_facs <- facs_results %>% 
  dplyr::rename(reporter = "type") %>% 
  filter(reporter != "noRE") %>% 
  mutate(dist = high_gate - low_gate) %>% 
  group_by(reporter) %>% 
  summarize(dist = geoMean(dist), mfi = geoMean(mfi))

plot_df <- left_join(mean_facs, n_df) %>% 
  group_by(reporter, dist, mfi) %>% 
  summarize(n = sum(n)) %>% 
  mutate(cluster = ifelse(reporter %in% c("RE57M", "RE57L"), "proximal",
                          ifelse(reporter %in% c("RE96", "RE61", "RE85"), "distal", "noRE-like")))

cor <- round(cor(plot_df$mfi, plot_df$n), 2)

plot_cor <- plot_df %>% 
  ggplot(aes(x = log10(mfi), y = n, color = factor(cluster, levels = c("proximal", "distal", "noRE-like")))) +
  geom_point(size = 0.5) +
  annotate(geom = "text", label = paste0("r=", cor), x = 3.3, y = 25, size = 6/2.8, color = "black") +
  scale_y_continuous(limits = c(0, NA), expand = c(0.1, 0.1)) +
  scale_x_continuous(expand = c(0.1, 0.1), breaks = c(2.5, 3, 3.5, 4)) +
  scale_color_manual(values = c("#90134D", "#F7A71C", "#676767"), name = "Cluster") 

fix <- set_panel_size(plot_cor, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9f_facs_interactions.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Plot per TF
tf_df <- mean_zfc %>% 
  filter(fdr <= 0.2) %>% 
  mutate(type = ifelse(mean_zfc <= 0, "rep", "act")) %>% 
  group_by(type, Gene) %>% 
  tally() %>% 
  mutate(Gene = str_remove(Gene, "_p[1-3]"))

order <- tf_df %>% 
  group_by(Gene) %>% 
  summarize(sum = sum(n)) %>% 
  arrange(-sum) %>% 
  select(Gene) %>% 
  unlist(use.names = FALSE)

tf_plot <- tf_df %>% 
  ggplot(aes(x = factor(Gene, levels = order), y = n, fill = type)) +
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values = c("#F7A71C", "#58BEBF")) +
  scale_y_continuous(name = "Number of interactions") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank())

fix <- set_panel_size(tf_plot, height = unit(2, "cm"), width = unit(14, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9d_tf_interactions.pdf", fix, dpi = 300,
       useDingbats=FALSE)
