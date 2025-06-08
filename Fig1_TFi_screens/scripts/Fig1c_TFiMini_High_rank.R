#Plots rank plot for TFiMini Library (High/Negative)
library(tidyverse)
library(egg)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  axis.ticks.x = element_blank(), axis.text.y = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(), 
                  strip.text = element_text(size = 6)))


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Designate controls
non_tf <- c("Dcp1a", "Rif1", "Brd4", "Lbr", "Stag2", "Nup62cl")

mle <- read.delim("./input_files/TFiMini_mle_HvN.gene_summary.txt") %>% 
  mutate(rank = rank(High.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = High.beta, pval = High.wald.p.value, fdr = High.wald.fdr, rank) %>% 
  mutate(identity = ifelse(gene %in% non_tf | is.na(prom), "cont", ifelse(pval >= 0.05, "ns", 
                                                                          ifelse(beta <= 0, "act", "rep"))))

#Calculate number of enriched/depleted TF genes
n <- mle %>% 
  filter(!gene %in% non_tf) %>% 
  filter(!is.na(prom)) %>% 
  group_by(identity) %>% 
  summarize(n = n())

#Plot results as rank plot
rank_plot <- mle %>% 
  ggplot(aes(x = rank, y = beta, color = identity)) +
  geom_point(size = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = mle[mle$rank <= 8,], aes(label = gene, color = identity), size = 6/2.8,
                  max.overlaps = 20) +
  geom_text_repel(data = mle[mle$rank >= 138,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text(data = n[n$identity == "act",], aes(x = 20, y = 1.4, label = paste0(n, "/119 TFs")),
            color = "#F7A71C", size = 6/2.8) +
  geom_text(data = n[n$identity == "rep",], aes(x = 120, y = -1.4, label = paste0(n, "/119 TFs")),
            color = "#58BEBF", size = 6/2.8) +
  labs(y = "Beta Score (High / Negative)") +
  scale_y_continuous(limits = c(-2.8, 2.8), breaks = c(-2, -1, 0, 1, 2)) +
  scale_color_manual(values = c("#F7A71C", "#000000", "#C6C6C6", "#58BEBF"))


fix <- set_panel_size(rank_plot, height = unit(4, "cm"), width = unit(4.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig1c_TFiMini_rank.pdf", fix, dpi = 300,
       useDingbats=FALSE)


