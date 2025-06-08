#Plots heatmaps comparing TFi High, TFi Low and TFiMini High target enrichment
library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)

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

#Specify non-tf controls
non_tf <- c("Xert", "Dcp1a", "Tsix", "Rif1", "Brd4", "Spen", "Dnmt1", "Trim28", "Kat8", "Eed", "Chd8",   
            "Msl1", "Msl2", "Kansl1", "Kansl3",  "Lbr", "Hnrnpu",  "Hnrnpk", "Bcorl1", "Mid2", "Map3k15", "Stag2",  
            "Jpx", "Bcor", "Eras", "Ripply1", "Tslrn1", "Nup62cl", "Xist", "Gm9785",  "Ftx",    
            "Nexmif", "Rlim")

#Read mageck mle output and mark activators/repressors
mle_Low <- read.delim("./input_files/TFi_mle_LvN.gene_summary.txt") %>% 
  mutate(rank = rank(Low.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = Low.beta, pval = Low.wald.p.value, fdr = Low.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep"))))


mle_High <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>% 
  mutate(rank = rank(High.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = High.beta, pval = High.wald.p.value, fdr = High.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep"))))

mle_TFiMini <- read.delim("./input_files/TFiMini_mle_HvN.gene_summary.txt") %>% 
  mutate(rank = rank(High.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = High.beta, pval = High.wald.p.value, fdr = High.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep"))))

TFimini_combi <- mle_TFiMini %>% 
  select(tu, beta_TFiMini = beta)

#Create lists of activators/repressors
Low_act <- unique(mle_Low[mle_Low$identity=="act",]$gene)
Low_rep <- unique(mle_Low[mle_Low$identity=="rep",]$gene)

High_act <- unique(mle_High[mle_High$identity=="act",]$gene)
High_rep <- unique(mle_High[mle_High$identity=="rep",]$gene)



#Create Heatmap
comb_high <- mle_High %>% 
  select(tu, gene, prom, beta_high = beta, identity_high = identity, rank_high = rank) 

comb_low <- mle_Low %>% 
  select(tu, gene, prom, beta_low = beta, identity_low = identity, rank_low = rank) 

#Merge all mle files together
#Structure by behavior in High and Low populations
combi <- left_join(comb_high, comb_low) %>% 
  filter(identity_high != "ns" | identity_low != "ns") %>% 
  mutate(identity_high = ifelse(identity_high != "ns" & beta_high <= 0, "highact", 
                                ifelse(identity_high != "ns" & beta_high >= 0, "highrep", "highns")),
         identity_low = ifelse(identity_low != "ns" & beta_low <= 0, "lowact", 
                               ifelse(identity_low != "ns" & beta_low >= 0, "lowrep", "lowns"))) %>%
  mutate(comb_ident = paste(identity_high, identity_low, sep = "_")) %>% 
  arrange(comb_ident, beta_high) %>% 
  left_join(TFimini_combi)

pheat_df <- combi %>% 
  select(tu, beta_high, beta_low, beta_TFiMini) %>% 
  pivot_longer(-tu, names_to = "comp", values_to = "beta") %>% 
  pivot_wider(names_from = tu, values_from = beta) %>% 
  column_to_rownames("comp")


annot_row <- combi %>% 
  select(tu, comb_ident) %>% 
  column_to_rownames("tu")

pheat_mat <- pheat_df %>% 
  as.matrix()

mat_col <- colorRampPalette(c("#F7A71C", "#FFFFFF", "#58BEBF"))(21)

rep_heat <- pheatmap(pheat_mat, 
                     annotation_col = annot_row,
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE,
                     border_color = "#666666",
                     fontsize = 6, 
                     gaps_col = c(19, 35, 40, 44, 62),
                     gaps_row = c(1, 2),
                     color = mat_col, 
                     breaks = seq(-1, 1.0952381, 0.0952381),
                     cellwidth = 0.22 / 0.0353, 
                     cellheight = 0.22 / 0.0353)

ggsave("./output_files/Fig1d_Screens_heat.pdf", rep_heat, dpi = 300,
       useDingbats=FALSE)      

