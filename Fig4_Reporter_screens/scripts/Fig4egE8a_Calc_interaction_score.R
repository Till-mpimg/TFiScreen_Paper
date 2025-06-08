#REPi analysis
library(tidyverse)
library(pheatmap)
library(egg)
library(gridExtra)
library(viridis)
library(ggrepel)

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


counts <- read.delim("./input_files/REPi.count_normalized.txt")
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")

#Excludes non-TF X-linked genes, as well as sgRNAs targeting REs
#Excludes Zfp866 because it almost completely overlaps with Gm20422
exclude <- c("RE57", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127", "FIREWACh", "Hnrnpk_p1", "Dcp1a_p1",
                      "RE46", "RE47", "RE49", "RE50", "Nup62cl_p1", "Lbr_p1", "RE51", "Zfp866_p2", "RE12", "RE53", 
             "Stag2_p3", "Rif1_p1", "Brd4_p1")

#Filters guides that have a low mean count across samples
filter_guides <-  counts %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("reporter", "frac", "replicate")) %>% 
  filter(frac == "Unsorted") %>% 
  group_by(sgRNA) %>% 
  summarize(mean_count = mean(counts)) %>% 
  filter(mean_count <= 150) %>% 
  select(sgRNA) %>% 
  unlist(use.names = FALSE)

#Calculates Log2Foldchange accross all guides and reporters
proc <- counts %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "counts") %>% 
  separate(sample, c("reporter", "frac", "replicate")) %>% 
  filter(frac != "Unsorted") %>% 
  filter(Gene != "NT") %>%
  filter(!sgRNA %in% filter_guides)


lfc <- proc %>% 
  pivot_wider(names_from = frac, values_from = counts) %>% 
  mutate(lfc = log2((High+0.01)/(Low+0.01)))

#Exports log2 fold changes to text file
write_delim(lfc, "./output_files/REPi_lfc.txt", delim = "\t")


#Calculate sig genes in noRE fraction
noRE_df <- lfc %>% 
  filter(reporter == "noRE") %>% 
  filter(!Gene %in% exclude) 

noRE_test <- expand.grid(Gene = unique(noRE_df$Gene), reporter = unique(noRE_df$reporter))

#Tests enrichment between High and Low fraction for each target with an unpaired t-test
noRE_test$pval <- mapply(function(a) {t.test(noRE_df[noRE_df$Gene==a,]$High, noRE_df[noRE_df$Gene==a,]$Low, 
                                                                                    var.equal = TRUE)$p.value},
                         noRE_test$Gene)
noRE_test$fdr <- p.adjust(noRE_test$pval) 

#Puts lfc and fdr in the same dataframe for plotting
plot_noRE_df <- left_join(noRE_df, noRE_test) %>% 
  group_by(Gene, fdr, pval) %>% 
  summarize(lfc = mean(lfc))

#Export noRE results to files
write_delim(plot_noRE_df, "./output_files/noRE_results.txt", delim = "\t")

#Plots volcano plot of targets affecting the noRE screens
noRE_plot <- plot_noRE_df %>% 
  ggplot(aes(x = lfc, y = -log10(fdr), group = Gene)) +
  geom_point(color = "#BCBCBC", size = 0.5) +
  geom_hline(aes(yintercept = -log10(0.2)), linetype = "dashed") +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(limits = c(-2.6, 2.6)) +
  geom_point(data = plot_noRE_df[plot_noRE_df$fdr <= 0.2,], size = 0.5) +
  geom_text_repel(data = plot_noRE_df[plot_noRE_df$fdr <= 0.2,], aes(label = Gene),
                  size = 6/2.8)


fix <- set_panel_size(noRE_plot, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8b_noRE_volcano.pdf", fix, dpi = 300,
       useDingbats=FALSE)

exclude_noRE <- noRE_test %>%  
  filter(fdr <= 0.2) %>% 
  select(Gene) %>% 
  unlist(use.names = FALSE)


#Plot PCA to compare reporters
lfc_wide <- lfc %>% 
  filter(!Gene %in% exclude) %>%
  filter(!Gene %in% exclude_noRE) %>% 
  unite("sample" ,reporter, replicate) %>% 
  group_by(Gene, sample) %>% 
  summarize(lfc = mean(lfc)) %>% 
  pivot_wider(names_from = sample, values_from = lfc) %>% 
  column_to_rownames("Gene")

pca <- prcomp(t(lfc_wide))
summary(pca)

pca_df <- data.frame(pca$x[,1:36]) %>%
  rownames_to_column("sample")
                      
#Writes principal components to text file
write_delim(pca_df,"./output_files/REPi_pca.txt" ,delim = "\t")

#Plots PC1/2 as scatter plot to show similarity between replicates
plot_pca <- data.frame(pca$x[,1:2]) %>% 
  rownames_to_column("sample") %>% 
  separate(sample, c("reporter", "replicate"), sep = "_") %>% 
  ggplot(aes(x = PC1, y = PC2, fill = factor(reporter, levels = c("noRE", level)))) +
  geom_point(shape = 21) +
  scale_x_continuous(name = "PC1 [30.0%]") +
  scale_y_continuous(name = "PC2 [18.7%]") +
  scale_fill_manual(values = c("#666666", viridis(n = 11, option = "turbo", begin = 0.1, end = 0.9))) +
  theme(legend.title = element_blank())


fix <- set_panel_size(plot_pca, height = unit(3, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4g_pca_lfc.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Sort reporter screens by kmeans
set.seed(22)
kmeans <- kmeans(t(lfc_wide), 3)

plot_pca_kmeans <- data.frame(pca$x[,1:2]) %>% 
  rownames_to_column("sample") %>% 
  left_join(data.frame(cluster = kmeans$cluster, sample = names(kmeans$cluster))) %>% 
  separate(sample, c("reporter", "replicate"), sep = "_") %>% 
  ggplot(aes(x = PC1, y = PC2, fill = as.factor(cluster))) +
  geom_point(shape = 21) +
  scale_fill_viridis(discrete = TRUE, option = "H", begin = 0.2, end = 0.8)  +
  theme(legend.title = element_blank())

fix <- set_panel_size(plot_pca_kmeans, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig4g_pca_kmeans.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Calculate normed lfc by subtracting noRE results
help_df <- lfc %>% 
  filter(reporter == "noRE") %>% 
  group_by(Gene) %>% 
  summarize(noRE = mean(lfc))

#Remove non-TF targets and significant interactors with the noRE line
norm_df <- lfc %>% 
  filter(!reporter == "noRE") %>% 
  filter(!Gene %in% exclude_noRE)%>% 
  filter(!Gene %in% exclude) %>% 
  select(-High, -Low) %>% 
  left_join(help_df) %>% 
  mutate(norm_lfc = lfc - noRE)


#Calculate interaction score (ZFC)
#Normalized log2 fold changes are zscore transformed and inversed (to get positive scores for positive interactions)
zfc <- norm_df  %>%
  group_by(reporter, replicate) %>% 
  mutate(zfc = as.numeric(-scale(norm_lfc))) %>% 
  ungroup()



#Calculate pval with one-sided t-test (Zscore scales around 0)
test_df <- expand.grid(Gene = unique(zfc$Gene), reporter = unique(zfc$reporter))

test_df$pval <- mapply(function(a,b) {t.test(zfc[zfc$Gene==a & zfc$reporter==b,]$zfc, mu = 0)$p.value},
                       test_df$Gene, test_df$reporter)

#Calculate FDR to correct for multiple testing
fdr_fun <- function(a) {
  fdr <- p.adjust(test_df[test_df$reporter==a,]$pval)
  Gene <- unique(test_df$Gene)
  
  rep_df <- data.frame(fdr, Gene) %>% 
    mutate(reporter = a)
  
  return(rep_df)
}


fdr_df <- bind_rows(lapply(unique(test_df$reporter), fdr_fun)) 

#Write to dataframe to output to text files
output_df <- left_join(zfc, test_df) %>% 
  left_join(fdr_df)

#Calculate mean results for later use
mean_df <- output_df %>% 
  group_by(Gene, reporter, pval, fdr) %>% 
  summarize(lfc = mean(lfc), norm_lfc = mean(norm_lfc), zfc = mean(zfc))

write_delim(output_df, "./output_files/REPi_results.txt",
            delim = "\t")

write_delim(mean_df, "./output_files/REPi_mean_results.txt",
            delim = "\t")

#Creates a matrix to use for pheatmap (only with significant genes)
filter_sig <- fdr_df %>% 
  filter(fdr <= 0.2) %>% 
  select(Gene) %>% 
  unique() %>% 
  unlist()

pheat_df <- zfc %>%
  filter(Gene %in% filter_sig) %>% 
  group_by(Gene, reporter) %>% 
  summarize(zfc = mean(zfc)) %>% 
  select(Gene, reporter, zfc) %>% 
  pivot_wider(names_from = reporter, values_from = zfc) %>% 
  mutate(Gene = str_remove(Gene, "_.*")) %>% 
  column_to_rownames("Gene") %>% 
  select(RE57L, RE57M, RE57R, RE58, RE61, RE85, RE93, RE95, RE96, RE97, RE127)

pheat_mat <- pheat_df %>% 
  as.matrix()

sig_df <- fdr_df %>% 
  filter(Gene %in% filter_sig) %>% 
  pivot_wider(names_from = reporter, values_from = fdr) %>% 
  mutate(Gene = str_remove(Gene, "_.*")) %>% 
  arrange(Gene) %>% 
  column_to_rownames("Gene") %>% 
  select(RE57L, RE57M, RE57R, RE58, RE61, RE85, RE93, RE95, RE96, RE97, RE127)

sig_mat <- sig_df %>% 
  as.matrix()

sig_mat[sig_mat <= 0.2] <- "*"
sig_mat[sig_mat != "*"] <- ""

mat_col <- colorRampPalette(c("#5BBFBF", "#FFFFFF", "#F7A81F"))(21)

rep_heat <- pheatmap(t(pheat_mat), 
                     border_color = "#C6C6C6",
                     cluster_rows = FALSE,
                     fontsize = 6, 
                     color = mat_col, 
                     cuttree_cols = 3,
                     gaps_row = c(4,10),
                     breaks = seq(-3, 3, 0.2857143),
                     display_numbers = t(sig_mat),
                     cellwidth = 0.2 / 0.0353, 
                     cellheight = 0.2 / 0.0353,
                     treeheight_col = 20, 
                     clustering_method = "ward.D",
                     legend = TRUE)

ggsave("./output_files/Fig4e_heatmap_zfc.pdf", rep_heat, dpi = 300,
       useDingbats=FALSE)      

