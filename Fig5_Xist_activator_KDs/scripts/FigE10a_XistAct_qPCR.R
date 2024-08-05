#Script analyzing qPCRs from the validation of the TFi Screen
library(tidyverse)
library(scales)
library(egg)
library(readxl)
library(EnvStats)
library(pheatmap)
library(gridExtra)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank(),
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Takes the standard .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC
raw <- list.files(path = "./input_files/", pattern = "TFiVal.*.xls", full.names = TRUE) %>% 
  lapply(read_xls, sheet = "Results", skip = 40) %>% 
  bind_rows %>%
  select("Sample Name", "Target Name", "CT") %>%
  dplyr::rename(Sample = "Sample Name", Gene = "Target Name") %>%
  filter(Sample != "NTC") %>%
  mutate(CT = ifelse(CT == "Undetermined", 40, CT)) %>% 
  transform(CT = as.numeric(CT)) %>%
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarize(CT = geoMean(CT)) %>%
  na.omit()

#Calculate relative expression to Arpo and Rrm2
analysis <- raw %>% 
  pivot_wider(names_from = Gene, values_from = CT) %>% 
  mutate(Cont = (Arpo + Rrm2) /2, Arpo = NULL, Rrm2 = NULL) %>%
  select(Sample, Cont, everything()) 

for( i in 3:length(analysis) ) {
  analysis[i] <- analysis[i] - analysis[2]
  analysis[i] <- 2^ - analysis[i]
  analysis[i] <- log2(analysis[i])
}

output <- analysis %>%
  select(-Cont) 

write.table(output, "./data/XistAct_qPCR_RelExp.txt",sep="\t",row.names = FALSE, quote = FALSE)


#Plots results - associate guide numbers to target genes
sgRNA <- c("LR15", "LR23", "TS14", "TS15", "TS16", "TS17", "TS18", "TS19", "TS20", "TS21", "TS22", "TS23", "TS24",
           "TS25")

Target <- c("NTC", "RE57", "Oct4", "Oct4", "Otx2", "Otx2", "Zic3", "Zic3", "Nfrkb", "Nfrkb", "Foxd3", "Foxd3", 
            "Nfe2l2","Nfe2l2") 
sg_num <- c("A", "A", rep(c("A", "B"), 6))


target_df <- data.frame(sgRNA, Target, sg_num)

res <- output %>% 
  separate(Sample, c("Rep", "sgRNA", "Med"), sep = "_") %>%
  pivot_longer(-c(1:3), names_to = "Gene", values_to = "Rel_Exp") %>% 
  left_join(target_df)

#Calculate pvalues for +dTAG/-dTAG comparisons
sgRNA <- unique(res$sgRNA)
test_df <- res %>% 
  select(Gene, sgRNA) %>% 
  unique()

test_df$pval <- mapply(function(a,b) {t.test(res[res$sgRNA==a & res$Gene==b & res$Med=="dTAG",]$Rel_Exp,
                                             res[res$sgRNA==a & res$Gene==b & res$Med=="DIFF",]$Rel_Exp,
                                             var.equal = TRUE)$p.value},
                       test_df$sgRNA, test_df$Gene)



#Plot log2 fold changes between dTAG and knockdown as heatmap
lfc_res <- res %>% 
  pivot_wider(names_from = Med, values_from = Rel_Exp) %>% 
  mutate(lfc = DIFF - dTAG)



#Plot the results as a heatmap
gene_level <- c("Xist", "Oct4", "Zic3", "Nfrkb", "Otx2", "Foxd3", "Nfe2l2", "Jpx", "Ftx", "Xert", "Rnf12", 
                "Linx", "Tsix", "Rex1", "Nanog", "Esrrb")

heat <- lfc_res %>% 
  group_by(sgRNA, Gene, Target, sg_num) %>% 
  summarize(lfc = mean(lfc)) %>% 
  pivot_wider(names_from = Gene, values_from = lfc)

pheat_mat <- heat  %>%
  arrange(factor(Target, levels = c("NTC", "RE57", "Oct4", "Zic3", "Nfrkb", "Otx2", "Foxd3", "Nfe2l2"))) %>%  
  unite(sample, Target, sg_num) %>% 
  ungroup() %>% 
  select(-sgRNA) %>% 
  column_to_rownames("sample") %>% 
  select(all_of(gene_level)) %>% 
  as.matrix()

pval_mat <- test_df %>% 
  left_join(heat[,1:3]) %>% 
  arrange(factor(Target, levels = c("NTC", "RE57", "Oct4", "Zic3", "Nfrkb", "Otx2", "Foxd3", "Nfe2l2"))) %>% 
  unite(sample, Target, sg_num) %>% 
  pivot_wider(names_from = Gene, values_from = pval) %>% 
  column_to_rownames("sample") %>%
  select(-sgRNA) %>% 
  select(all_of(gene_level)) %>% 
  as.matrix()

pval_mat[pval_mat <= 0.05] <- "*"
pval_mat[pval_mat != "*"] <- ""

mat_col <- colorRampPalette(c("#F7A71C", "#FFFFFF", "#58BEBF"))(21)

target_heat <- pheatmap(pheat_mat, 
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE,
                        gaps_row = c(1,2, 4, 6, 8, 10,12),
                        gaps_col = c(1, 7, 11, 13),
                        border_color = "#666666",
                        fontsize = 6, 
                        color = mat_col, 
                        display_numbers = pval_mat,
                        breaks = seq(-3, 3.2857143, 0.2857143), 
                        cellwidth = 0.25 / 0.0353, 
                        cellheight = 0.25 / 0.0353)


ggsave("./output_files/FigE10a_qpcr.pdf", target_heat, dpi = 300,
       useDingbats=FALSE)

