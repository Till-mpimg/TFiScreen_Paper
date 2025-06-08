#Compares mle data from TFi and TFiMini
library(tidyverse)
library(egg)
library(gridExtra)
library(ggvenn)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title = element_text(size = 6),
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

non_tf <- c("Xert", "Dcp1a", "Tsix",  "Spen", "Dnmt1", "Kat8", "Eed", "Chd8",   
            "Msl1", "Msl2", "Kansl1", "Kansl3",  "Lbr", "Hnrnpu",  "Hnrnpk", "Bcorl1", "Mid2", "Map3k15", "Stag2",  
            "Jpx", "Bcor", "Eras", "Ripply1", "Tslrn1", "Nup62cl", "Xist", "Gm9785",  "Ftx",    
            "Nexmif", "Rlim")

#Load mle for High/Negative comparisons from both screens
mle_TFi <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, TFi.beta = High.beta, TFi.pval = High.wald.p.value, TFi.fdr = High.wald.fdr)


mle_TFiMini <- read.delim("./input_files/TFiMini_mle_HvN.gene_summary.txt") %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, TFiMini.beta = High.beta, TFiMini.pval = High.wald.p.value, TFiMini.fdr = High.wald.fdr)

mle <- inner_join(mle_TFi, mle_TFiMini)

#Put both mle results in the same dataframe
sig_df <- mle %>% 
  mutate(sig = ifelse(TFi.pval <= 0.05 & TFiMini.pval <= 0.05, "both", 
                      ifelse(TFi.pval <= 0.05, "TFi",
                             ifelse(TFiMini.pval <= 0.05, "TFiMini", "ns")))) %>% 
  filter(!gene %in% non_tf)

#Export comparison between the screens as .txt file
write_delim(sig_df, "./output_files/TFiLib_comp.txt", delim= "\t")

pearson <- round(cor(mle$TFi.beta, mle$TFiMini.beta), 2)

#Plot scatter plot to compare betascores
scatter_plot <- sig_df %>% 
  ggplot(aes(x = TFi.beta, y = TFiMini.beta, color = sig)) +
  geom_point(size = 0.5) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1), expand = c(0, 0.1)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1), expand = c(0, 0.1)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  scale_color_manual(values = c("#EA6A08", "#E6E022", "#F7A71C", "#C6C6C6")) +
  geom_text(aes(x = -0.5, y = 0.5, label = paste0("r=", pearson)), color = "black", size = 6/2.8)

fix <- set_panel_size(scatter_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3h_TFi_comp_scatter.pdf", fix, dpi = 300,
       useDingbats=FALSE)



#Create lists of activators/repressors
TFi_act <- unique(sig_df[sig_df$TFi.pval<=0.05 & sig_df$TFi.beta<=0,]$gene)
TFi_rep <- unique(sig_df[sig_df$TFi.pval<=0.05 & sig_df$TFi.beta>=0,]$gene)

TFiMini_act <- unique(sig_df[sig_df$TFiMini.pval<=0.05 & sig_df$TFiMini.beta<=0,]$gene)
TFiMini_rep <- unique(sig_df[sig_df$TFiMini.pval<=0.05 & sig_df$TFiMini.beta>=0,]$gene)


#Create Venn Diagrams to visualize overlapping significant genes
venn_list <- list("TFi_act" = TFi_act, "TFi_rep" = TFi_rep, "TFiMini_act" =  TFiMini_act, "TFiMini_rep" =  TFiMini_rep)

act_venn <- ggvenn(venn_list, c("TFi_act", "TFiMini_act"), show_percentage = FALSE, auto_scale = TRUE, set_name_size = 6/2.8,
                   text_size = 6/2.8, stroke_size = 0.5)

fix <- set_panel_size(act_venn, height = unit(1, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3i_TFi_comp_act_venn.pdf", fix, dpi = 300,
       useDingbats=FALSE)

rep_venn <- ggvenn(venn_list, c("TFi_rep", "TFiMini_rep"), show_percentage = FALSE, auto_scale = TRUE, set_name_size = 6/2.8,
                   text_size = 6/2.8, stroke_size = 0.5)

fix <- set_panel_size(rep_venn, height = unit(1, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE3i_TFi_comp_rep_venn.pdf", fix, dpi = 300,
       useDingbats=FALSE)
