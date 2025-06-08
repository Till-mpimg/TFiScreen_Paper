#Do analysis on the reporter screen output to compare to TFi results
library(tidyverse)
library(egg)
library(gridExtra)

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
mle_Low_TFi <- read.delim("./input_files/TFi_mle_LvN.gene_summary.txt")
mle_High_TFi <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt")
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")


#Compares all interactions to mle
mean_zfc <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(mean_zfc = mean(zfc))

sig_zfc <- mean_zfc %>%  
  left_join(mle_Low_TFi) %>% 
  left_join(mle_High_TFi) %>%
  select(Gene, reporter, mean_zfc, Low.beta, High.beta)

#Do a cumulative distribution ranked by the different beta scores
low_order <- sig_zfc %>% 
  ungroup() %>% 
  select(Gene, Low.beta) %>% 
  unique() %>% 
  arrange(Low.beta) %>% 
  mutate(low_order = row_number())

high_order <- sig_zfc %>% 
  ungroup() %>% 
  select(Gene, High.beta) %>% 
  unique() %>% 
  arrange(High.beta) %>% 
  mutate(high_order = row_number())



cum_df <- sig_zfc %>% 
  left_join(low_order) %>% 
  left_join(high_order) %>% 
  pivot_longer(c(6:7), names_to = "comp", values_to = "rank") %>% 
  group_by(comp, reporter) %>% 
  arrange(rank) %>% 
  mutate(cumul = cumsum(mean_zfc)) %>% 
  select(comp, rank, cumul, reporter)

sup_tab <- cum_df %>% 
  pivot_wider(names_from = comp, values_from = cumul)

write_delim(sup_tab, "./output_files/SuppTable_cumulScore.txt", delim = "\t")


origin_point <- expand.grid(comp = unique(cum_df$comp), rank = 0, cumul = 0, 
                            reporter = unique(cum_df$reporter))


plot_df <- cum_df %>% 
  rbind(origin_point)

##CREATE PERMUTATED DATASETS
scramble_df <- sig_zfc %>% 
  select(Gene, reporter, mean_zfc)
scramble_gene <- unique(scramble_df$Gene)

n <- 10000
# List to store permuted datasets
permuted_datasets <- list()

# Perform permutations
set.seed(123)  # Set seed for reproducibility
for (i in 1:n) {
  scramble_order <- data.frame(Gene = sample(scramble_gene)) %>% 
    mutate(rank = row_number()) 
  
  loop_df <- scramble_df %>% 
    left_join(scramble_order) %>% 
    group_by(reporter) %>% 
    arrange(rank) %>% 
    mutate(cumul = cumsum(mean_zfc)) %>% 
    select(rank, cumul, reporter)
  
  permuted_datasets[[i]] <- loop_df
}

#Check if means are significantly different from permuted datasets
actual_means <- cum_df %>% 
  group_by(comp, reporter) %>% 
  summarize(mean_score = mean(cumul))

permuted_means <- lapply(permuted_datasets, 
                         function(df) {
                           out <- df %>% 
                             group_by(reporter) %>% 
                             summarize(mean_score = mean(cumul))
                          return(out)
                         })
means <- bind_rows(permuted_means)

actual_means$pval <- mapply(function(a, b) {
  fil_means <- means %>% 
    filter(reporter==a)
  p_value <- mean(means[means$reporter==a,]$mean_score >= 
                    actual_means[actual_means$comp==b & actual_means$reporter==a,]$mean_score)
  return(p_value)
}, actual_means$reporter, actual_means$comp)

actual_means$mean_delta <- mapply(function(a, b) {
  fil_means <- means %>% 
    filter(reporter==a)
  mean_delta <- mean(actual_means[actual_means$comp==b & actual_means$reporter==a,]$mean_score - means$mean_score)
  return(mean_delta)
}, actual_means$reporter, actual_means$comp)
  
ci_df <- bind_rows(permuted_datasets) %>% 
  group_by(rank, reporter) %>% 
  summarize(perc_highest = quantile(cumul, 0.99), perc_high = quantile(cumul, 0.95))

cum_plot <- plot_df %>% 
  ggplot(aes(x = rank)) +
  facet_wrap(~factor(reporter, levels = level)) + 
  geom_hline(aes(yintercept = 0), color = "black") +
  geom_step(data = ci_df, aes(y = perc_highest), color = "#676767", linetype = "dotted", alpha = 0.5) +
  geom_step(data = ci_df, aes(y = perc_high), color = "#ABABAB", linetype = "dotted", alpha = 0.5) +
  geom_step(aes(y = cumul, color = comp)) +
  scale_color_manual(values = c("#B21A29", "#F28A41")) +
  scale_y_continuous(name = "Cumulative interaction score [REPi screen]", breaks = c(0, 8, 16)) +
  scale_x_continuous(name = "Gene rank [TFi screen]", breaks = c(0, 50, 100)) +
  geom_text(data = actual_means, aes(label = pval, x = 100, y = 12, color = comp), size = 4/2.8)

fix <- set_panel_size(cum_plot, height = unit(1, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4h_cumulative_interactionscore.pdf", fix, dpi = 300,
       useDingbats=FALSE)
