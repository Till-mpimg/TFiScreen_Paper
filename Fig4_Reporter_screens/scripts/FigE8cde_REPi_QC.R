#This script plots QC for REPi Screens
library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), 
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 6),  
                  strip.text = element_text(size = 6), strip.background = element_blank(),
                  legend.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Read counts from screen fractions
counts <- read.delim("./output_files/REPi.count_normalized.txt")
level <- c("noRE", "RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")


#Calculates distribution width of the fractions
width_df <- counts %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low) %>% 
  separate(sample, c("reporter", "frac", "replicate"))


width_plot <- width_df %>%
  ggplot(aes(x = width, y = factor(reporter, levels = rev(level)), color = frac)) +
  facet_wrap(~replicate) +
  geom_point(size = 1) +
  labs(x = "Distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1.3, 3.4) +
  scale_color_manual(values = c("#F7A71C", "#58BEBF", "#666666"))

fix <- set_panel_size(width_plot, height = unit(3, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8c_REPi_log2_width.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Compute only for NTCs
ntc_width_df <- counts %>%
  filter(Gene == "NT") %>%
  pivot_longer(c(-sgRNA, -Gene), names_to = "sample", values_to = "counts") %>%
  group_by(sample) %>%
  summarize(width_low = log2(quantile(counts, 0.1)), width_high = log2(quantile(counts, 0.9))) %>%
  mutate(width = width_high - width_low)  %>% 
  separate(sample, c("reporter", "frac", "replicate"))



ntc_width_plot <- ntc_width_df %>%
  ggplot(aes(x = width, y = factor(reporter, levels = rev(level)), color = frac)) +
  facet_wrap(~replicate) +
  geom_point(size = 1) +
  labs(x = "Distribution width (log2)") +
  theme(axis.title.y = element_blank()) +
  xlim(1.3, 3.4) +
  scale_color_manual(values = c("#F7A71C", "#58BEBF", "#666666"))

fix <- set_panel_size(ntc_width_plot, height = unit(3, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE8d_REPi_NTC_width.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Pearson between all samples
filter_guides <-  counts %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "counts") %>%
  separate(sample, c("reporter", "frac", "replicate")) %>% 
  filter(frac == "Unsorted") %>% 
  group_by(sgRNA) %>% 
  summarize(mean_count = mean(counts)) %>% 
  filter(mean_count <= 150) %>% 
  select(sgRNA) %>% 
  unlist(use.names = FALSE)

pearson_df <- counts %>% 
  filter(Gene != "NT") %>% 
  select(-Gene) %>% 
  filter(!sgRNA %in% filter_guides) %>% 
  column_to_rownames("sgRNA")

cor_mat <- cor(pearson_df)


#Produce pearson heatmaps in a loop
loop_vec <- c("noRE", "RE57L", "RE57M", "RE57R",
              "RE58", "RE61", "RE85", "RE93", "RE95",
              "RE96", "RE97", "RE127")
pearson_fun <- function(a) {
  rows_to_keep <- grepl(a, rownames(cor_mat))
  cols_to_keep <- grepl(a, colnames(cor_mat))
  
  filtered_matrix <- cor_mat[rows_to_keep, cols_to_keep]
  rownames(filtered_matrix) <- sub(paste0("_", a), "", rownames(filtered_matrix))
  colnames(filtered_matrix) <- sub(paste0("_", a), "", colnames(filtered_matrix))
  
  mat_col <- colorRampPalette(c("blue", "white", "red"))(11)
  
  cor_heat <- pheatmap(filtered_matrix, 
                       border_color = "#000000",
                       fontsize = 6, 
                       color = mat_col, 
                       breaks = seq(0, 1, 0.1),
                       cellwidth = 0.22 / 0.0353, 
                       cellheight = 0.22 / 0.0353,
                       treeheight_row = 0, 
                       treeheight_col = 0, 
                       legend = TRUE)
  
  ggsave(paste0("./output_files/FigE8e_REPi_", a,  "_pearson.pdf"),cor_heat, dpi = 300,
         useDingbats=FALSE)
}

sapply(loop_vec, pearson_fun)

