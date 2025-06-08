#Plots rank plot for TFiMini control populations
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


#Read mle data for PlasmidLib vs Selected comp
mle_Sel <- read.delim("./input_files/TFiMini_mle_SvP.gene_summary.txt") %>% 
  mutate(rank = rank(Selection.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = Selected.beta, pval = Selected.wald.p.value, fdr = Selected.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(beta <= 0, "depleted", "enriched")))

#Count number of enriched targets

n <- mle_Sel %>% 
  group_by(tu) %>% 
  select(tu, identity) %>% 
  ungroup() %>% 
  group_by(identity) %>% 
  summarize(n = n())

#Plot effects of targets on cell growth during selection as rank plot
rank_plot <- mle_Sel %>% 
  ggplot(aes(x = rank, y = beta, color = identity)) +
  geom_point(size = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = mle_Sel[mle_Sel$rank <= 3,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text_repel(data = mle_Sel[mle_Sel$rank >= 137 & mle_Sel$identity =="enriched",], 
                  aes(label = gene, color = identity), size = 6/2.8) +
  geom_text(data = n[n$identity == "depleted",], aes(x = 40, y = -0.8, label = paste0(n, "/142 Targets")),
            color = "#1C6CF7", size = 6/2.8) +
  geom_text(data = n[n$identity == "enriched",], aes(x = 100, y = 0.8, label = paste0(n, "/142 Targets")),
            color = "#DE2D26", size = 6/2.8) +
  labs(y = "Beta Score (Selected/Plasmid Library)") +
  ylim(-1.6, 1.6) +
  scale_color_manual(values = c("#1C6CF7", "#DE2D26", "#C6C6C6"))


fix <- set_panel_size(rank_plot, height = unit(3, "cm"), width = unit(4.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3f_TFiMini_Selected_rank.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Read mle data for 2iSL vs Selected comp
mle_2iSL <- read.delim("./input_files/TFiMini_mle_Sv2iSL.gene_summary.txt") %>% 
  mutate(rank = rank(X2iSL.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = X2iSL.beta, pval = X2iSL.wald.p.value, fdr = X2iSL.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(beta <= 0, "depleted", "enriched")))

#Count number of enriched targets
n <- mle_2iSL %>% 
  group_by(tu) %>% 
  select(tu, identity) %>% 
  ungroup() %>% 
  group_by(identity) %>% 
  summarize(n = n())

#Plot effects of targets on cell growth during tissue culture as rank plot
rank_plot <- mle_2iSL %>% 
  ggplot(aes(x = rank, y = beta, color = identity)) +
  geom_point(size = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = mle_2i[mle_2i$rank <= 1,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text(data = n[n$identity == "depleted",], aes(x = 20, y = -0.8, label = paste0(n, "/142 Targets")),
            color = "#1C6CF7", size = 6/2.8) +
  geom_text(data = n[n$identity == "enriched",], aes(x = 120, y = 0.8, label = paste0("0/142 Targets")),
            color = "#DE2D26", size = 6/2.8) +
  labs(y = "Beta Score (2iSL/Selected)") +
  ylim(-1.6, 1.6) +
  scale_color_manual(values = c("#1C6CF7", "#C6C6C6" ))


fix <- set_panel_size(rank_plot, height = unit(3, "cm"), width = unit(4.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3f_TFiMini_2iSL_rank.pdf", fix, dpi = 300,
       useDingbats=FALSE)



