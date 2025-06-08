#Plots rank plot for TFi Library (Low/Negative)
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

#Designate controls (+High confidence controls)
non_tf <- c("Xert", "Dcp1a", "Tsix", "Rif1", "Brd4", "Spen", "Dnmt1", "Trim28", "Kat8", "Eed", "Chd8",   
            "Msl1", "Msl2", "Kansl1", "Kansl3",  "Lbr", "Hnrnpu",  "Hnrnpk", "Bcorl1", "Mid2", "Map3k15", "Stag2",  
            "Jpx", "Bcor", "Eras", "Ripply1", "Tslrn1", "Nup62cl", "Xist", "Gm9785",  "Ftx",    
            "Nexmif", "Rlim")
high_conf <- c("Xert", "Tsix", "Rif1", "Spen", "Brd4", "Msl1", "Msl2", "Kat8", "Jpx", "Ftx", "Rlim")

#Load Low/Negative mle input
mle <- read.delim("./input_files/TFi_mle_LvN.gene_summary.txt") %>% 
  mutate(rank = rank(Low.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = Low.beta, pval = Low.wald.p.value, fdr = Low.wald.fdr, rank) %>% 
  mutate(identity = ifelse(gene %in% high_conf, "cont",
                           ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep")))))

#Calculate number of enriched/depleted tf genes
n <- mle %>% 
  group_by(gene) %>% 
  slice_min(pval) %>% 
  filter(!gene %in% non_tf) %>% 
  select(gene, identity) %>% 
  ungroup() %>% 
  group_by(identity) %>% 
  summarize(n = n())


  
#Plots results as a rank plot
#High confidence controls + controls with p<0.05 are colored in dark gray
rank_plot <- mle %>% 
  ggplot(aes(x = rank, y = beta, color = identity)) +
  geom_point(size = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.2) +
  geom_text_repel(data = mle[mle$rank <= 6,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text_repel(data = mle[mle$rank >= 908,], aes(label = gene, color = identity), size = 6/2.8) +
  geom_text(data = n[n$identity == "act",], aes(x = 200, y = 0.8, label = paste0(n, "/570 TFs")),
            color = "#F7A71C", size = 6/2.8) +
  geom_text(data = n[n$identity == "rep",], aes(x = 700, y = -0.8, label = paste0(n, "/570 TFs")),
            color = "#58BEBF", size = 6/2.8) +
  labs(y = "Beta Score (Low/Negative)") +
  ylim(-1.6, 1.6) +
  scale_color_manual(values = c("#F7A71C", "#000000", "#C6C6C6", "#58BEBF"))


fix <- set_panel_size(rank_plot, height = unit(3, "cm"), width = unit(4.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE2h_TFi_rank_Low.pdf", fix, dpi = 300,
       useDingbats=FALSE)


