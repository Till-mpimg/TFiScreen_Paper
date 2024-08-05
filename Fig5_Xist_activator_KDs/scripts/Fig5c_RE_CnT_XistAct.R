#Xist activator CnT analysis
library(tidyverse)
library(Rsubread)
library(gridExtra)
library(EnvStats)
library(egg)
library(pheatmap)



theme_set(theme_classic() +
            theme(legend.title = element_text(size = 6), plot.title = element_text(hjust = 0.5),
                  axis.title = element_text(size = 6), panel.border = element_rect(size = 0.5, color = "black", fill = NA),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  legend.text = element_text(size = 6),
                  strip.background = element_blank(), strip.text = element_text(size = 6),
                  legend.key.size = unit(0.5,"line")))


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


bam_dir <- "./data/"

setwd(bam_dir)
temp = list.files(pattern=paste0(".*_dedup.bam$"), full.names = TRUE)
candidate_re <- "./input_files/re_total.saf"

feature_counts <- featureCounts(temp, annot.ext = candidate_re, isPairedEnd = TRUE,
                                allowMultiOverlap = TRUE)


counts_df <- data.frame(feature_counts$counts)
names(counts_df) <- gsub(x = names(counts_df), pattern = "_dedup.bam", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "\\.", replacement = "_")
names(counts_df) <- gsub(x = names(counts_df), pattern = "XX_D2_", replacement = "")


total_reads <- feature_counts$stat %>%
  dplyr::select(-Status) %>%
  summarize_all(sum) %>%
  mutate_all(funs(. / 1000000))

CPM_df <- sweep(counts_df, 2, as.vector(t(unlist(total_reads))), FUN = '/') %>%
  rownames_to_column("RE") %>% 
  pivot_longer(-RE, names_to = "sample", values_to = "cpm") %>%
  separate(sample, c("mark", "line", "rep"), sep = "_")

supp_table <- CPM_df %>% 
  pivot_wider(names_from = rep, values_from = cpm)

write_delim(supp_table, "./output_files/CnT_Xist_REs_CPM.txt", delim = "\t")

sig_df <- expand.grid(RE = unique(CPM_df$RE), mark = unique(CPM_df$mark), line = c("sgZic3", "sgNfrkb", "sgOtx2", "sgFoxd3"))

sig_fun <- function(a, b, c) {
  t.test(CPM_df[CPM_df$RE==a & CPM_df$mark==b & CPM_df$line=="sgNT",]$cpm, 
         CPM_df[CPM_df$RE==a & CPM_df$mark==b & CPM_df$line==c,]$cpm, var.equal = TRUE)$p.value
}

sig_df$pval <- mapply(sig_fun, sig_df$RE, sig_df$mark, sig_df$line)

avg_counts <- counts_df  %>%
  rownames_to_column("RE") %>% 
  pivot_longer(-RE, names_to = "sample", values_to = "count") %>%
  separate(sample, c("mark", "rep"), sep = "_") %>% 
  group_by(RE, mark) %>% 
  summarize(avg_counts = mean(count), min_counts = min(count), max_counts = max(count)) %>% 
  filter(avg_counts <= 5) %>% 
  mutate(filt_name = paste0(RE, "-", mark))
  

lfc_df <- CPM_df %>% 
  group_by(RE, mark, line) %>% 
  summarize(cpm = mean(cpm)) %>% 
  pivot_wider(names_from = line, values_from = cpm) %>% 
  mutate(across(where(is.numeric), ~ (. + 0.01)/(sgNT + 0.01))) %>% 
  select(-sgNT) %>% 
  pivot_longer(-c(RE, mark), names_to = "line", values_to = "fc") %>% 
  mutate(lfc = log2(fc), filt_name = paste0(RE, "-", mark))

sup_table_lfc <- lfc_df %>% 
  filter(RE %in% c("RE_57L", "RE_57M", "RE_57R", "RE_58", "RE_61", "RE_85", 
                   "RE_93", "RE_95", "RE_96", "RE_97", "RE_127")) %>% 
  left_join(sig_df) %>% 
  mutate(lfc = ifelse(filt_name %in% avg_counts$filt_name, NA, lfc), 
         pval = ifelse(filt_name %in% avg_counts$filt_name, NA, pval)) %>% 
  select(RE, mark, line, lfc, pval)

write_delim(sup_table_lfc, "/output_files/CnT_Xist_REs_LFC.txt", delim = "\t")


pheat_df <- lfc_df %>%
  mutate(name = paste0(mark, "_", line), lfc = ifelse(filt_name %in% avg_counts$filt_name, NA, lfc)) %>%
  ungroup() %>% 
  select(RE, name, lfc) %>%
  pivot_wider(names_from = name, values_from = lfc)  %>% 
  filter(RE %in% c("RE_57L", "RE_57M", "RE_57R", "RE_58", "RE_61", "RE_85", 
                   "RE_93", "RE_95", "RE_96", "RE_97", "RE_127")) %>% 
  column_to_rownames("RE")

pval_df <- sig_df %>% 
  filter(RE %in% c("RE_57L", "RE_57M", "RE_57R", "RE_58", "RE_61", "RE_85", 
                   "RE_93", "RE_95", "RE_96", "RE_97", "RE_127")) %>% 
  mutate(name = paste0(mark, "_", line), filt_name = paste0(RE, "-", mark)) %>% 
  mutate(pval = ifelse(filt_name %in% avg_counts$filt_name, 1, pval)) %>% 
  select(RE, name, pval) %>%
  pivot_wider(names_from = name, values_from = pval)%>% 
  column_to_rownames("RE")

pval_mat <- pval_df %>% 
  as.matrix()

pval_mat[pval_mat <= 0.05] <- "*"
pval_mat[pval_mat != "*"] <- ""


pheat_mat <- pheat_df %>% 
  as.matrix()

pheat_mat <-  pheat_mat[match(c("RE_57L", "RE_57M", "RE_57R", "RE_58", "RE_61", "RE_85", 
                                "RE_93", "RE_95", "RE_96", "RE_97", "RE_127"), rownames(pheat_mat)), ]


pval_mat <-  pval_mat[match(c("RE_57L", "RE_57M", "RE_57R", "RE_58", "RE_61", "RE_85", 
                              "RE_93", "RE_95", "RE_96", "RE_97", "RE_127"), rownames(pval_mat)), ]

pval_mat <-  pval_mat[,match(colnames(pheat_mat), colnames(pval_mat))]

mat_col <- colorRampPalette(c("#F7A81F", "#FFFFFF", "#5BBFBF"))(21)

heat <- pheatmap(t(pheat_mat), 
                     border_color = "#C6C6C6",
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize = 6, 
                     color = mat_col, 
                     gaps_row = c(4, 8),
                     breaks = seq(-2.5, 2.5, 0.2380952),
                     display_numbers = t(pval_mat),
                     cellwidth = 0.2 / 0.0353, 
                     cellheight = 0.2 / 0.0353,
                     legend = TRUE)

ggsave("./output_files/Fig5c_CnT_heat.pdf", heat, useDingbats = FALSE)


