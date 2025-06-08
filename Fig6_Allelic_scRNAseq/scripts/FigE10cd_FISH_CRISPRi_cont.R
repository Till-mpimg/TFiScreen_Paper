#Plots FISH output as Fig. 5e-g and Extended Fig. E10d
library(tidyverse)
library(EnvStats)
library(egg)
library(gridExtra)
library(viridis)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]

setwd(wd)


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6), legend.position = "none",axis.title.x = element_blank(), 
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))



#Load segmentation output
fish_total <- read.delim("./output_files/CRISPRi_cont_analysis_out.txt") %>% 
  mutate(sample = str_remove(sample, "SP107_")) %>% 
  separate(sample, c("rep", "line", "day"), sep = "_") %>% 
  mutate(comp = ifelse(line %in% c("AMC7", "LR15"),"SP107", "DEL"))


fish_sum <- read.delim("./output_files/CRISPRi_cont_analysis_sum.txt") %>% 
  mutate(sample = str_remove(sample, "SP107_")) %>% 
  separate(sample, c("rep", "line", "day"), sep = "_") %>% 
  mutate(comp = ifelse(line %in% c("AMC7", "LR15"),"SP107", "DEL"))


#Plots % of Xist-cells
test_df <- data.frame(comp = unique(fish_sum$comp)) %>% 
  mutate(cont = ifelse(comp == "DEL", "A3", "LR15"), line = ifelse(comp == "DEL", "F9", "AMC7"))

test_df$pval <- mapply(function(a, b) {t.test(fish_sum[fish_sum$line==a,]$xist_perc,
                                           fish_sum[fish_sum$line==b,]$xist_perc, 
                                             var.equal = TRUE, paired = TRUE)$p.value},
                       test_df$cont, test_df$line)

#Plot percentage of Xist positive cells
xist_plot <- fish_sum %>%
  ggplot(aes(x = factor(line, levels = c("A3", "LR15", "F9", "AMC7")), y = xist_perc)) +
  facet_wrap(~comp, scales = "free_x") +
  stat_summary(aes(fill = line), 
                   geom = "bar", fun = "mean", position = position_dodge(width = 1)) +
  stat_summary(aes(group = line, y = bi_perc), 
               geom = "bar", fun = "mean", position = position_dodge(width = 1), color = "black", fill = NA) +
  geom_point(aes(fill = line),
             size = 0.7, position = position_dodge(width = 1), shape = 21) +
  scale_y_continuous(limits = c(0, 1), name = "%Xist-positive cells") +
  geom_text(data = na.omit(test_df[test_df$pval<=0.05,]), aes(label = "*", y = 0.85, x = 1.5), 
            size = 8/2.8) +
  scale_fill_manual(values = c("#CDCDCD", "#AA0023", "#AA0023", "#CDCDCD"))

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE10c_CRISPRi_cont_FISH_perc.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Plots amount of signal per cloud (bg corrected)
cloud_df <- fish_total %>% 
  filter(n_xist %in% c(1,2))


#Plot Xist level per cell downsampled to match the lowest replicate
min_cloud_df <- cloud_df %>% 
  group_by(line, rep) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(line) %>% 
  summarize(min = min(n))

cloud_downsample <- cloud_df %>% 
  left_join(min_cloud_df) %>% 
  group_by(line, rep) %>% 
  sample_n(size = first(min), replace = FALSE) %>% 
  ungroup()


test_df$p_clouds <- mapply(function(a,b) {wilcox.test(cloud_downsample[cloud_downsample$line==a,]$cy5_norm_sub,
                                                       cloud_downsample[cloud_downsample$line==b,]$cy5_norm_sub)$p.value},
                           test_df$cont, test_df$line)


downsample_plot <- cloud_downsample %>%
  ggplot(aes(x = factor(line, levels = c("A3", "LR15", "F9", "AMC7")), y = log10(cy5_norm_sub))) +
  facet_wrap(~comp, scales = "free_x") +
  geom_dotplot(aes(fill = line), 
               binaxis = "y", stackdir = "center", stackratio = 0.8, 
               binwidth = 0.05, dotsize = 1.5, stroke = 0.1, 
               position = position_dodge(width = 0.8)) +
  stat_summary(aes(group = line), geom = "crossbar", fun = "geoMean", color = "black", lwd = 0.25, width = 0.25,
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(name = "Xist intensity/Cloud [a.u. - log10]", limits = c(1.5, 5), breaks = c(2,3,4)) +
  geom_text(data = na.omit(test_df[test_df$p_clouds<=0.05,]), aes(label = "*", y = 4.8),
            size = 8/2.8)+
  scale_fill_manual(values = c("#CDCDCD", "#AA0023", "#AA0023", "#CDCDCD"))

fix <- set_panel_size(downsample_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE10d_CRISPRi_cont_FISH_intensity.pdf", fix, dpi = 300,
       useDingbats=FALSE)

