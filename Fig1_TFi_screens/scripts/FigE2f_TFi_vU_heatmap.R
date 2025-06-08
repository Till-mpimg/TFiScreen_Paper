#Plots log2foldchanges for individual replicates vs the Unsorted fraction
library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)
library(EnvStats)

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

#Reads HvN mle file and selects significant hits
mle <- read.delim("./input_files/TFi_mle_HvN.gene_summary.txt") %>% 
  filter(Gene != "NT") %>% 
  mutate(rank = rank(High.beta, ties.method = "first")) %>% 
  separate(Gene, c("gene", "prom"), remove = FALSE, sep = "_") %>% 
  select(tu = Gene, gene, prom, beta = High.beta, pval = High.wald.p.value, fdr = High.wald.fdr, rank) %>% 
  mutate(identity = ifelse(pval >= 0.05, "ns", ifelse(gene %in% non_tf, "cont", ifelse(beta <= 0, "act", "rep"))))

act <- mle$tu[mle$identity=="act"]
rep <- mle$tu[mle$identity=="rep"] 

#Sets order according to rank in HvN comparison
order <- mle %>% 
  filter(identity %in% c("act", "rep")) %>% 
  mutate(order_col = ifelse(identity == "act", rank, 100000/rank)) %>% 
  arrange(order_col) %>% 
  select(tu) %>% 
  unlist(use.names = FALSE)

#Reads in counts for all significant hits
counts <- read.delim("./input_files/TFi.count_corr_normalized.txt") %>%
  filter(Gene %in% c(act, rep))

#Calculates fold change per guide compared to the Unsorted fraction
fc_df <- counts %>% 
  pivot_longer(-c(1:2), names_to = "frac", values_to = "count") %>%
  filter(frac != "PlasmidLib") %>% 
  separate(frac, c("rep", "frac"), sep = "_") %>% 
  filter(frac %in% c("Unsorted", "High", "Negative")) %>% 
  pivot_wider(names_from = frac, values_from = count) %>% 
  mutate(High = High/Unsorted, Negative = Negative/Unsorted)

#Summarizes fold change per promoter and converts to log2
sum_df <- fc_df %>% 
  group_by(rep, Gene) %>% 
  summarize(High = geoMean(High), Negative = geoMean(Negative)) %>% 
  pivot_longer(-c(1:2), names_to = "frac", values_to = "fc") %>% 
  mutate(lfc = log2(fc)) %>% 
  mutate(identity = ifelse(Gene %in% act, "act", "rep"))

#Creates a matrix to use for pheatmap
pheat_df <- sum_df %>% 
  unite(sample, rep, frac) %>% 
  select(-fc) %>% 
  pivot_wider(names_from = sample, values_from = lfc) %>% 
  arrange(factor(Gene, levels = order)) %>% 
  column_to_rownames("Gene") %>% 
  select(R1_High, R2_High, R1_Negative, R2_Negative, identity)

pheat_idents <- data.frame(idents = pheat_df$identity)
rownames(pheat_idents) <- rownames(pheat_df)
  
pheat_mat <- pheat_df %>% 
  select(-identity) %>% 
  as.matrix()

mat_col <- colorRampPalette(c("#F7A71C", "#FFFFFF", "#58BEBF"))(21)

#Plots mean LFCs per replicate
rep_heat <- pheatmap(pheat_mat, 
         annotation_row = pheat_idents,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         border_color = "#666666",
         fontsize = 6, 
         gaps_col = 2,
         gaps_row = 26,
         color = mat_col, 
         breaks = seq(-1, 1.0952381, 0.0952381),
         cellwidth = 0.22 / 0.0353, 
         cellheight = 0.22 / 0.0353)

ggsave("./output_files/FigE2f_TFi_vU_heat.pdf", rep_heat, dpi = 300,
       useDingbats=FALSE)      
