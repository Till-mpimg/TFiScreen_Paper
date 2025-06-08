#Plots FISH output as Fig. 5e-g and Extended Fig. E10d
library(tidyverse)
library(EnvStats)
library(egg)


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

fish_total <- read.delim("./output_files/sgAct_FISH_analysis_out.txt")
  ))
fish_sum <- read.delim("./output_files/sgAct_FISH_analysis_sum.txt")

#Plots % of Xist-cells
test_df <- data.frame(guide = unique(fish_sum$guide))

test_df$pval <- mapply(function(a) {t.test(fish_sum[fish_sum$guide=="sgNT",]$xist_perc,
                                             fish_sum[fish_sum$guide==a,]$xist_perc, 
                                             var.equal = TRUE, paired = TRUE)$p.value},
                       test_df$guide)

target_levels = c("sgNT", "sgZic3", "sgNfrkb", "sgOtx2", "sgFoxd3")

#Plot percentage of Xist positive cells
xist_plot <- fish_sum %>%
  ggplot(aes(x = factor(guide, levels = target_levels), y = xist_perc)) +
  stat_summary(aes(fill = factor(guide, levels = target_levels)), 
                   geom = "bar", fun = "mean", position = position_dodge(width = 1)) +
  stat_summary(aes(group = factor(guide, levels = target_levels), y = bi_perc), 
               geom = "bar", fun = "mean", position = position_dodge(width = 1), color = "black", fill = NA) +
  geom_point(aes(fill = factor(guide, levels = target_levels)),
             size = 0.7, position = position_dodge(width = 1), shape = 21) +
  scale_y_continuous(limits = c(0, 1), name = "%Xist-positive cells") +
  scale_x_discrete(name = "Guide") +
  geom_text(data = na.omit(test_df[test_df$pval<=0.05,]), aes(label = "*", y = 0.85), 
            size = 8/2.8)

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig5e_sgAct_FISH_perc.pdf", fix, dpi = 300,
       useDingbats=FALSE)



#Plots amount of signal per cloud (bg corrected)
cloud_df <- fish_total %>% 
  filter(n_xist %in% c(1,2))


#Plot Xist level per cell downsampled to match the lowest replicate
min_cloud_df <- cloud_df %>% 
  group_by(guide, rep) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(guide) %>% 
  summarize(min = min(n))

cloud_downsample <- cloud_df %>% 
  left_join(min_cloud_df) %>% 
  group_by(guide, rep) %>% 
  sample_n(size = dplyr::first(min), replace = FALSE) %>% 
  ungroup()


test_df$p_clouds <- mapply(function(a) {wilcox.test(cloud_downsample[cloud_downsample$guide=="sgNT",]$cy5_norm_sub,
                                                       cloud_downsample[cloud_downsample$guide==a,]$cy5_norm_sub)$p.value},
                               test_df$guide)


downsample_plot <- cloud_downsample %>% 
  ggplot(aes(x = factor(guide, levels = target_levels), 
             y = log10(cy5_norm_sub))) +
  geom_dotplot(aes(fill = factor(guide, target_levels)), 
               binaxis = "y", stackdir = "center", stackratio = 0.8, 
               binwidth = 0.05, dotsize = 1.5, stroke = 0.1, 
               position = position_dodge(width = 0.8)) +
  stat_summary(aes(group = guide), geom = "crossbar", fun = "geoMean", color = "black", lwd = 0.25, width = 0.5,
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(name = "Xist intensity/Cloud [a.u. - log10]", limits = c(1, 5), breaks = c(2,3,4)) +
  scale_x_discrete(name = "Guide")  +
  geom_text(data = na.omit(test_df[test_df$p_clouds<=0.05,]), aes(label = "*", y = 4.7),
            size = 8/2.8)

fix <- set_panel_size(downsample_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig5f_sgAct_FISH_signal_intensity_downsample.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#plot intensity by replicate for supplement
test_df_rep <- expand.grid(guide = unique(fish_sum$guide), rep = unique(fish_sum$rep))



test_df_rep$signal_p <- mapply(function(a, b) {wilcox.test(cloud_df[cloud_df$guide=="sgNT" & cloud_df$rep==b,]$cy5_norm_sub,
                                                    cloud_df[cloud_df$guide==a & cloud_df$rep==b,]$cy5_norm_sub)$p.value},
                               test_df_rep$guide, test_df_rep$rep)

rep_plot <- cloud_df %>% 
  ggplot(aes(x = factor(guide, levels = target_levels), 
             y = log10(cy5_norm_sub))) +
  facet_wrap(~rep) +
  geom_dotplot(aes(fill = factor(guide, levels = target_levels)), 
               binaxis = "y", stackdir = "center", stackratio = 0.8, 
               binwidth = 0.05, dotsize = 1.5, stroke = 0.1, 
               position = position_dodge(width = 0.8)) +
  stat_summary(aes(group = guide), geom = "crossbar", fun = "geoMean", color = "black", lwd = 0.25, width = 0.5,
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(name = "Xist intensity/Cloud [a.u. - log10]", limits = c(1, 5), breaks = c(2,3,4)) +
  scale_x_discrete(name = "Guide")  +
  geom_text(data = na.omit(test_df_rep[test_df_rep$signal_p<=0.05,]), aes(label = "*", y = 4.7),
            size = 8/2.8)

fix <- set_panel_size(rep_plot, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9i_sgAct_FISH_signal_intensity_reps.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Plot in correlation to total Xist KD
mean_signal <- cloud_df %>%
  group_by(rep, guide) %>% 
  summarize(xist_signal = geoMean(cy5_norm_sub))

nt_signal <- mean_signal %>% 
  filter(guide == "sgNT") %>% 
  left_join(fish_sum[,c(1,2,7)]) %>% 
  select(rep, nt_signal = xist_signal, nt_perc = xist_perc)

xist_kd <- mean_signal %>% 
  left_join(fish_sum[,c(1,2,7)]) %>% 
  left_join(nt_signal) %>% 
  filter(guide != "sgNT") %>% 
  mutate(norm_signal = log10(xist_signal)- log10(nt_signal) , norm_perc = xist_perc - nt_perc)


xist_norm <- xist_kd %>%
  ggplot(aes(x = norm_signal, y = norm_perc, 
             color = factor(guide, levels = c("sgZic3", "sgNfrkb", "sgOtx2", "sgFoxd3")))) +
  geom_point(aes(shape = factor(guide, levels = c("sgZic3", "sgNfrkb", "sgOtx2", "sgFoxd3"))), size = 1.5) +
  scale_y_continuous(limits = c(-0.6, 0.05)) +
  scale_x_continuous(limits = c(-1, 0.05)) +
  scale_shape_manual(values = c(16, 15, 16, 15)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") 

fix <- set_panel_size(xist_norm, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig5g_sgAct_FISH_signal_comp.pdf", fix, dpi = 300,
       useDingbats=FALSE)

