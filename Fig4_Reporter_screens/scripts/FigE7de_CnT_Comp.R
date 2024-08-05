library(Rsubread)
library(tidyverse)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Read RE locations and FACS data
re_total <- read.delim("./input_files/re_total.bed",
                       col.names = c("chr", "start", "end", "RE", "trash1", "trash2", "trash3", "trash4", "trash5"),
                       header = FALSE) %>% 
  select(RE, chr, start, end)
  
facs_timecourse <- read.delim("./output_files/Firewach_Timecourse.txt")


mfi_df <- facs_timecourse %>% 
  group_by(rep, day, type, run) %>% 
  dplyr::summarize(mfi = geoMean(GFP)) 

help_df <- mfi_df %>% 
  filter(type %in% c("noRE", "neg")) %>% 
  pivot_wider(names_from = type, values_from = mfi) %>% 
  mutate(noRE = noRE - neg)

lfc_df <- mfi_df  %>% 
  filter(!type %in% c("noRE", "neg"))  %>%
  left_join(help_df) %>%  
  mutate(mfi = mfi - neg) %>% 
  mutate(lfc = log2(mfi/noRE)) %>% 
  select(-noRE, -neg)

mean_mfi <- lfc_df %>% 
  filter(day %in% c("D0", "D2", "D4")) %>% 
  group_by(day, type) %>% 
  summarize(lfc = mean(lfc))

#Read BAM files and count reads in REs
re_saf <- re_total %>% 
  transmute(GeneID = RE, Chr = chr, Start = start, End = end, Strand = ".")

bam_dir <- ("./input_files/Gjaltema_bam/")
bam_list <- list.files(bam_dir, pattern = ".bam$", full.names = TRUE)

feature_counts <- featureCounts(bam_list, annot.ext = re_saf, isPairedEnd = TRUE)

counts_df <- data.frame(feature_counts$counts)
names(counts_df) <- gsub(x = names(counts_df), pattern = "\\.", replacement = "_")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_merged_bam", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "_TX1072_XXdXic", replacement = "")
names(counts_df) <- gsub(x = names(counts_df), pattern = "CUTnTag_", replacement = "")


#Calculate CPM
total_reads <- feature_counts$stat %>%
  dplyr::select(-Status) %>%
  summarize_all(sum) %>%
  mutate_all(funs(. / 1000000))

CPM_df <- sweep(counts_df, 2, as.vector(t(unlist(total_reads))), FUN = '/') %>%
  rownames_to_column("type") %>% 
  pivot_longer(-type, names_to = "sample", values_to = "cpm") %>% 
  separate(sample, c("mark", "day"), sep = "_") %>% 
  mutate(day = toupper(day), type = str_remove(type, "_"))

#Combines the sequencing and facs data in a single dataframe
plot_df <- left_join(CPM_df, mean_mfi) %>% 
  mutate(cluster = ifelse(type %in% c("RE57L", "RE57M", "RE57R", "RE58"), "proximal",
                          ifelse(type %in% c("RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127"), "distal",
                          ifelse(type %in% c("RE12", "RE46", "RE47", "RE50", "RE51", "RE52", "RE53"), "repressive", "other")))

cor_df <- expand.grid(mark = c("ATAC", "H3K27ac"), day = c("D0", "D2", "D4"))

#Calculates correlation between the two modalities
cor_fun <- function(a, b) {
  fun_df <- plot_df[plot_df$mark == a & plot_df$day == b,]
  cor <- cor(log2(fun_df$cpm + 1), fun_df$lfc,method = "spearman" )
}
cor_df$cor <- mapply(cor_fun, cor_df$mark, cor_df$day)

scatter_comp <- plot_df %>% 
  ggplot(aes(y = log2(cpm+1), x = lfc, color = factor(cluster, levels = c("proximal", "distal", "repressive", "other")))) +
  facet_wrap(mark~day) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_point(size = 0.5) +
  geom_text(data = cor_df, aes(label = paste0("r=", round(cor, 1)), x = -2, y = 4), size = 5/2.8, color = "black") +
  scale_x_continuous(limits = c(-3.2, 5), name = "Rel. MFI to empty Reporter [log2]") +
  scale_color_manual(values = c("#90134D", "#F7A71C", "5BBFBF", "#666666"), name = "Cluster") +
  scale_y_continuous(name = "Signal strength [log2(CPM + 1)]")

fix <- set_panel_size(scatter_comp, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE7de_Cnt_comp.pdf", fix, dpi = 300,
       useDingbats=FALSE)

