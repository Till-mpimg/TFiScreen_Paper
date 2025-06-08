library(tidyverse)
library(egg)
library(gridExtra)
library(drc)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]
setwd(wd)

xist_bin_col <- colorRampPalette(c("#ABABAB", "#93134D"))(10)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

cell_df <- read.delim("./output_files/scXist_filtered_cell_table.txt")

#Bin cells according to Xist expression
bins_reps <- cell_df %>% 
  dplyr::select(cell, Xist, rep, sample) %>% 
  unique() %>% 
  group_by(rep, sample) %>% 
  mutate(xist_bin = factor(ntile(Xist, n = 10)))

bins_reps %>% group_by(rep, sample, xist_bin) %>% tally()

bin_reps_df <- cell_df %>% 
  left_join(bins_reps)





#Check if alleles are different
#Plot Xist expression and allelic ratio per bin
sum_bin_reps_df <- bin_reps_df %>%
  filter(gene != "Xist") %>% 
  group_by(Xist, cell, rep, sample, xist_bin) %>% 
  summarize(sum_xi = sum(count_xi), sum_xa = sum(count_xa)) %>% 
  mutate(chrX_ratio = sum_xi / (sum_xi + sum_xa))

#Normalize Xist expression to bin 10
norm_factor_df <- sum_bin_reps_df %>% 
  filter(xist_bin == 10) %>% 
  mutate(Xist =  exp(Xist) - 1) %>% 
  group_by(xist_bin, sample, rep) %>% 
  summarize(norm_factor = median(Xist)) %>% 
  ungroup() %>% 
  dplyr::select(-xist_bin)

xist_norm_df <- sum_bin_reps_df %>% 
  mutate(Xist =  exp(Xist) - 1) %>% 
  group_by(xist_bin, sample, rep) %>% 
  summarize(Xist = median(Xist)) %>% 
  left_join(norm_factor_df) %>%  
  mutate(norm_Xist = Xist / norm_factor)

plot_bin_df <- sum_bin_reps_df %>% 
  group_by(xist_bin, sample, rep) %>% 
  summarize(chrX_ratio = median(chrX_ratio)) %>% 
  left_join(xist_norm_df)


#Do model
split_df <- plot_bin_df %>%
  group_by(sample, rep) %>%
  group_split()

sample_names <- plot_bin_df %>%
  group_by(sample, rep) %>%
  group_keys() %>%
  pull(sample, rep)

# Fit model per sample
grouped_models <- split_df %>%
  set_names(sample_names) %>%
  map(function(df) {
    tryCatch({
      model <- drm(chrX_ratio ~ norm_Xist, 
                   data = df, 
                   fct = LL.4(fixed = c(NA, 0, 0.5, NA)),
                   start = c(1, 0.25))
      df$predicted <- predict(model)
      df
    }, error = function(e) {
      df$predicted <- NA
      df
    })
  })

# Combine results into one dataframe
plot_bin_pred_df <- bind_rows(grouped_models)


chrX_normXist <- plot_bin_pred_df %>% 
  ggplot() +
  facet_wrap(sample~rep) +
  geom_line(aes(x = norm_Xist, y = predicted)) +
  geom_point(aes(x = norm_Xist, y = chrX_ratio), size = 1, color = "#93134D") +
  scale_x_continuous(limits = c(0, 1), name = "norm. Xist expression [to max bin]") +
  scale_y_continuous(limits = c(0, 0.55), name = "Allelic Ratio [chrX]", breaks = c(0, 0.25, 0.5))

model_fix <- set_panel_size(chrX_normXist, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(model_fix)
ggsave("./output_files/FigE10p_XistTitration_full_chrX_doseresponse_reps.pdf", model_fix, 
       useDingbats=FALSE)
