#Plots log2foldchanges for individual guides
library(tidyverse)
library(egg)
library(gridExtra)
library(EnvStats)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.position = "none",
                  axis.text.y = element_text(size = 6), axis.title.x = element_blank(), 
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

non_tf <- c("Xert", "Dcp1a", "Tsix", "Rif1", "Brd4", "Spen", "Dnmt1", "Trim28", "Kat8", "Eed", "Chd8",   
            "Msl1", "Msl2", "Kansl1", "Kansl3",  "Lbr", "Hnrnpu",  "Hnrnpk", "Bcorl1", "Mid2", "Map3k15", "Stag2",  
            "Jpx", "Bcor", "Eras", "Ripply1", "Tslrn1", "Nup62cl", "Xist", "Gm9785",  "Ftx",    
            "Nexmif", "Rlim")

#Reads HvN mle file and selects significant hits
mle <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>% 
  filter(Gene != "NT") %>% 
  mutate(rank = rank(High.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = High.beta, pval = High.wald.p.value, fdr = High.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep"))))

order <- mle %>% 
  filter(identity != "ns") %>% 
  arrange(rank) %>% 
  select(tu) %>% 
  unlist(use.names = FALSE)

#Reads in normalized counts
counts <- read.delim("./input_files/TFi.count_corr_normalized.txt") 

#Calculates fold change per guide compared to the Unsorted fraction
fc_df <- counts %>% 
  pivot_longer(-c(1:2), names_to = "frac", values_to = "count") %>%
  filter(frac != "PlasmidLib") %>% 
  separate(frac, c("rep", "frac"), sep = "_") %>% 
  filter(frac %in% c("High", "Negative")) %>% 
  pivot_wider(names_from = frac, values_from = count) %>% 
  mutate(fc = High/Negative)

#Summarizes replicates together and filters for significant hits
sum_df <- fc_df %>% 
  group_by(sgRNA, Gene) %>% 
  summarize(fc = geoMean(fc)) %>% 
  filter(Gene %in% mle[mle$identity!="ns",]$tu) %>%
  mutate(lfc = log2(fc), identity = ifelse(Gene %in% mle[mle$beta<=0,]$tu, "act", "rep")) 

#Plots fold change of individual guides for activators and repressors
guide_plot_act <- sum_df %>% 
  filter(identity == "act") %>% 
  separate(Gene, c("gene_id", "prom"), remove = FALSE) %>% 
  mutate(cont = ifelse(gene_id %in% non_tf, "cont", "tf")) %>% 
  ggplot(aes(x = factor(Gene, levels = order), y = lfc, color = cont, fill = cont)) +
  geom_point(size = 0.3) +
  geom_boxplot(alpha = 0.5, color = "black", size = 0.1, outlier.shape = NA, coef = FALSE,
               width = 0.7) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(limits = c(-3.1, 3.1), name = "Fold Change High/Negative [Log2]") +
  scale_color_manual(values = c("#666666", "#F7A71C")) +
  scale_fill_manual(values = c("#666666", "#F7A71C")) 

fix <- set_panel_size(guide_plot_act, height = unit(5, "cm"), width = unit(7, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE2g_TFi_Guides_act.pdf", fix, dpi = 300,
       useDingbats=FALSE)

guide_plot_rep <- sum_df %>% 
  filter(identity == "rep") %>%
  separate(Gene, c("gene_id", "prom"), remove = FALSE) %>% 
  mutate(cont = ifelse(gene_id %in% non_tf, "cont", "tf")) %>% 
  ggplot(aes(x = factor(Gene, levels = order), y = lfc, color = cont, fill = cont)) +
  geom_point(size = 0.3) +
  geom_boxplot(alpha = 0.5, color = "black", size = 0.1, outlier.shape = NA, coef = FALSE,
               width = 0.7) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_y_continuous(limits = c(-3.1, 3.1), name = "Fold Change High/Negative [Log2]") +
  scale_color_manual(values = c("#666666", "#58BEBF")) +
  scale_fill_manual(values = c("#666666", "#58BEBF")) 

fix <- set_panel_size(guide_plot_rep, height = unit(5, "cm"), width = unit(4.2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE2g_TFi_Guides_rep.pdf", fix, dpi = 300,
       useDingbats=FALSE)
