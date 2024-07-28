#Plots rank plot for TFi Library (Selected/PlasmidLibrary)
library(tidyverse)
library(egg)
library(gridExtra)
library(ggrepel)

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

#Read mle target enrichment and detect significantly enriched or depleted targets
mle <- read.delim("./input_files/TFi_mle_SvP.gene_summary.txt") %>% 
  mutate(rank = rank(Selected.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = Selected.beta, pval = Selected.wald.p.value, fdr = Selected.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(beta <= 0, "depleted", "enriched")))

#Calculate number of enriched or depleted targets
n <- mle %>% 
  group_by(tu) %>% 
  select(tu, identity) %>% 
  ungroup() %>% 
  group_by(identity) %>% 
  summarize(n = n())

#Plot effects of targets on cell growth during selection as rank plot
rank_plot <- mle %>% 
  ggplot(aes(x = rank, y = beta, color = identity)) +
  geom_point(size = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = mle[mle$rank <= 3,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text_repel(data = mle[mle$rank >= 909,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text(data = n[n$identity == "enriched",], aes(x = 700, y = -0.8, label = paste0(n, "/911 Promoters")),
            color = "#DE2D26", size = 6/2.8) +
  geom_text(data = n[n$identity == "depleted",], aes(x = 200, y = 0.8, label = paste0(n, "/911 Promoters")),
            color = "#1C6CF7", size = 6/2.8) +
  labs(y = "Beta Score (Selected/Plasmid Library)") +
  ylim(-1, 1) +
  scale_color_manual(values = c("#1C6CF7", "#DE2D26", "#C6C6C6"))


fix <- set_panel_size(rank_plot, height = unit(3, "cm"), width = unit(4.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1j_TFi_Selection_rank.pdf", fix, dpi = 300,
       useDingbats=FALSE)


