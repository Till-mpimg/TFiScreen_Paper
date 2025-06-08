#Do analysis on the reporter screen output
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
mle_TFi <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt")
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")


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

sig_zfc <- mean_zfc %>%  
  left_join(mle_TFi) %>% 
  mutate(type = ifelse(Gene %in% act, "act",
                       ifelse(Gene %in% rep, "rep", "ns")))


#Produce scatterplot
cor_fun = function(x){ 
  a <- sig_zfc[sig_zfc$reporter==x,]
  cor.p <- cor(a[a$fdr <= 0.2,]$mean_zfc, a[a$fdr <= 0.2,]$High.beta, method = "pearson")
  return(cor.p)
}

cor_df <- data.frame(reporter = unique(sig_zfc$reporter))
cor_df$cor <- sapply(cor_df$reporter, cor_fun)

scatter_zfc <- sig_zfc %>% 
  filter(fdr <= 0.2) %>% 
  ggplot(aes(x = mean_zfc, y = High.beta, color = type)) +
  facet_wrap(~factor(reporter, levels = level)) +
  geom_point(size = 0.45) +
  scale_x_continuous(limits = c(-4.85, 4.85), breaks = c(-4, -2, 0, 2, 4), name = "Enrichment score [REPi screen]") +
  scale_y_continuous(limits = c(-1, 1), name = "Beta score [TFi screen]") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_color_manual(values = c("#F7A71C", "#666666", "#58BEBF")) +
  geom_text(data = cor_df, 
            aes(x = 2.75, y = 0.75, label = paste0("r= ", round(cor, 2))), color = "black", size = 6/2.8)

fix <- set_panel_size(scatter_zfc, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8g_REPi_TFi_correlation.pdf", fix, dpi = 300,
       useDingbats=FALSE)

