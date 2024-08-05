#Compare reporter screen output depending on motif occurrences
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
motifs <- read.delim("./output_files/motifs_re_clean.txt")

#Calculate mean ZFC per interaction
mean_zfc <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(mean_zfc = mean(zfc))

#Filter tfs with any significant interactions
sig_vec <- mean_zfc %>% 
  separate(Gene, c("gene", "prom")) %>% 
  filter(fdr <= 0.2) %>% 
  ungroup() %>% 
  select(gene) %>% 
  unique() %>% 
  unlist(use.names = FALSE)

#Filter RE/TF combinations with a motif
fil_motifs <- motifs %>% 
  filter(RE %in% zfc$reporter) %>% 
  mutate(motif_tfs = str_extract(name, "\\(.*\\)")) %>% 
  mutate(motif_tfs = str_remove_all(motif_tfs, "[\\(\\)]")) %>% 
  mutate(motif_tfs = str_split(pattern = ";", motif_tfs)) %>% 
  unnest(motif_tfs) %>% 
  mutate(gene = str_to_title(motif_tfs), reporter = RE) %>% 
  select(reporter, gene, max_score)

#Filter enrichment scores for tfs with a motif in any of the tested reporters
#Assign repressors and activators from TFi Screen
#Calculate score as absolute value
motif_zfc <- mean_zfc %>% 
  separate(Gene, c("gene", "prom")) %>% 
  filter(gene %in% fil_motifs$gene) %>% 
  full_join(fil_motifs) %>% 
  filter(gene %in% sig_vec) %>% 
  mutate(motif = ifelse(is.na(max_score), FALSE, TRUE), 
         type = ifelse(gene %in% sig_TFi[sig_TFi$TFi.beta<=0 & sig_TFi$sig!="ns",]$gene , "act", 
                  ifelse(gene %in% sig_TFi[sig_TFi$TFi.beta>=0 & sig_TFi$sig!="ns",]$gene, "rep", "ns"))) %>% 
  select(-max_score) %>% 
  unique() %>% 
  mutate(abs_zfc = abs(mean_zfc)) %>% 
  mutate(dir = ifelse(mean_zfc > 0, "act", "rep"))

#Calculate number of factors per fraction
motif_zfc %>% ungroup() %>%  select(gene, type, reporter, motif) %>% unique() %>% group_by(type, motif) %>%  tally()

#Calculate pvalue for activators and repressors (and n.s.) seperately
pval_df <- data.frame(type = c("act", "rep", "ns"))

pval_df$pval <- sapply(pval_df$type, function(x) {
  motif_true <- motif_zfc[motif_zfc$motif == TRUE & motif_zfc$type == x,]$abs_zfc
  motif_false <- motif_zfc[motif_zfc$motif == FALSE & motif_zfc$type == x,]$abs_zfc
  wilcox.test(motif_true, motif_false)$p.value
})

#Plot as cumulative distribution separated by results of the TFi screen
ecdf_plot <- motif_zfc %>%
  ggplot(aes(x = abs_zfc)) +
  facet_wrap(~type) +
  stat_ecdf(size = 0.25, aes(color = motif)) +
  geom_text(data = pval_df, 
            aes(y = 0.5, x = 3, label = paste0("p=", round(pval, digits = 2))), size = 6/2.8, color = "black") +
  scale_color_manual(values = c("#666666", "#F7A71C")) +
  scale_x_continuous(name = "Enrichment score [absolute]", breaks = c(0, 2, 4)) +
  scale_y_continuous(name = "Cumulative distribution", breaks = c(0, 0.5, 1))

fix <- set_panel_size(ecdf_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9e_reporter_ecdf_motif_sep.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Plot as cumulative distribution for all the factors together
motif_true <- motif_zfc[motif_zfc$motif == TRUE,]$abs_zfc
motif_false <- motif_zfc[motif_zfc$motif == FALSE,]$abs_zfc
p_all <- wilcox.test(motif_true, motif_false)$p.value

ecdf_plot <- motif_zfc %>%
  ggplot(aes(x = abs_zfc)) +
  stat_ecdf(size = 0.25, aes(color = motif)) +
  annotate("text", 
            y = 0.5, x = 3, label = paste0("p=", round(p_all, digits = 2)), size = 6/2.8, color = "black") +
  scale_color_manual(values = c("#666666", "#F7A71C")) +
  scale_x_continuous(name = "Enrichment score [absolute]", breaks = c(0, 2, 4)) +
  scale_y_continuous(name = "Cumulative distribution", breaks = c(0, 0.5, 1))

fix <- set_panel_size(ecdf_plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9e_reporter_ecdf_motif_all.pdf", fix, dpi = 300,
       useDingbats=FALSE)



