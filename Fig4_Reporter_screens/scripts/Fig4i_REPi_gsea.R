#Do analysis on the reporter screen output
library(tidyverse)
library(egg)
library(gridExtra)
library(fgsea)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)
theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

zfc <- read.delim("./input_files/REPi_results.txt") 
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")

#Compares all interactions to mle
mean_zfc <- zfc %>% 
  group_by(Gene, reporter, fdr) %>% 
  summarize(mean_zfc = mean(zfc)) %>% 
  mutate(Gene = str_extract(Gene, "[A-Za-z0-9]*"))

sig_genes <- mean_zfc %>% 
  filter(fdr <= 0.2) %>% 
  ungroup() %>% 
  select(Gene) %>% 
  unique() %>% 
  unlist()

#Select tf groups
anova <- read.delim("./input_files/ANOVA_sex_tfs.txt") %>% 
  filter(gene %in% sig_genes)

cluster <- read.delim("./input_files/Timecourse_TF_clusters.txt") %>% 
  filter(gene %in% sig_genes)

naive <- cluster[cluster$cluster=="naive",]$gene
transient1 <- cluster[cluster$cluster=="form2",]$gene 
transient2 <- cluster[cluster$cluster=="form1",]$gene
transient3 <- cluster[cluster$cluster=="prime2",]$gene  
committed <- cluster[cluster$cluster=="prime1",]$gene 
xx_biased <- anova[anova$bias=="XX",]$gene

#Set list for GSEA analysis
clus_list <- list(naive = naive,
                  transient1 = transient1,
                  transient2 = transient2,
                  transient3 = transient3,
                  committed = committed,
                  xx_biased = xx_biased
)

#Put genes in matrix for analysis
zfc_mat <- mean_zfc %>% 
  filter(Gene %in% sig_genes) %>% 
  select(-fdr) %>% 
  pivot_wider(names_from = reporter, values_from = mean_zfc) %>% 
  column_to_rownames("Gene")

gsea_fun <- function(x) {
  rep_vec <- zfc_mat[,x]
  names(rep_vec) <- row.names(zfc_mat)
  
  gsea_res <- fgsea(
    pathways = clus_list,
    stats = rep_vec, 
    minSize = 5, 
    maxSize = 500,
    nperm = 1000,
    scoreType = "pos"
  )

  out_df <- data.frame(cluster = gsea_res$pathway, pval = gsea_res$pval, padj = gsea_res$padj, NES = gsea_res$NES,
                       ES = gsea_res$NES, reporter = x)
  return(out_df)
  
}

#Perform GSEA for every reporter and put results in a single dataframe
#Bin p values for plotting
gsea_out <- bind_rows(lapply(level, gsea_fun)) %>% 
  mutate(pval_bin = ifelse(pval <= 0.01, "0.01", 
                           ifelse(pval <= 0.1, "0.1", 
                                  ifelse(pval <= 0.2, "0.2",
                                  ifelse(pval <= 0.5, "0.5", ">0.5")))))
mat_col <- colorRampPalette(c("#FFFFFF", "#fcbba1", "#a50f15"))(10)


#Plot figure 4i
dotplot <- gsea_out %>% 
  ggplot(aes(x = factor(cluster, levels = c("naive", "transient1", "transient2", "transient3", "committed", "xx_biased")),
             y = factor(reporter, levels = rev(level)), fill = NES, size = pval_bin)) +
         scale_fill_stepsn(colors = mat_col, n.breaks = 11, limits = c(0, 3)) +
         scale_size_manual(values = c(">0.5" = 0.5, "0.5" = 0.875, "0.2" = 1.25, "0.1" = 1.625, "0.01" = 2)) +
         geom_point(shape = 21)

fix <- set_panel_size(dotplot, height = unit(3, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4i_REPi_gsea.pdf", fix, dpi = 300,
       useDingbats=FALSE)

