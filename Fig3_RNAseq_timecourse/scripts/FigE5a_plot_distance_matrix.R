library(tidyverse)
library(viridis)
library(pheatmap)


theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Read distance matrix
dist_mat <- read.delim("./output_files/distance_matrix.txt")

#Separates by sex and plots with pheatmap
pheat_df <- dist_mat %>% 
  select(1:5) %>% 
  filter(!str_detect(X, "pbulk")) %>% 
  separate(X, c("timepoint", "sex"))

pheat_XX <- pheat_df %>% 
  filter(sex == "XX") %>% 
  column_to_rownames("timepoint") %>% 
  select(clust_E3.5_pbulk, clust_E4.5_pbulk, clust_E5.5_pbulk, clust_E6.5_pbulk) %>% 
  as.matrix()

pheat_XO <- pheat_df %>% 
  filter(sex == "XO") %>% 
  column_to_rownames("timepoint") %>% 
  select(clust_E3.5_pbulk, clust_E4.5_pbulk, clust_E5.5_pbulk, clust_E6.5_pbulk) %>% 
  as.matrix()

mat_col <- colorRampPalette(c("#C42A19", "#FFFFFF", "#4F73B7"))(9)

XX_heat <- pheatmap(pheat_XX, 
                     border_color = "#FFFFFF",
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize = 6, 
                     color = mat_col, 
                     display_numbers = TRUE,
                     breaks = seq(40, 130, 10),
                     cellwidth = 1 / 0.0353, 
                     cellheight = 0.5 / 0.0353,
                     legend = TRUE)

ggsave("./output_files/FigE5a_XX_dist_heat.pdf", XX_heat, dpi = 300,
       useDingbats=FALSE)      

XO_heat <- pheatmap(pheat_XO, 
                    border_color = "#FFFFFF",
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    fontsize = 6, 
                    color = mat_col, 
                    display_numbers = TRUE,
                    breaks = seq(40, 130, 10),
                    cellwidth = 1 / 0.0353, 
                    cellheight = 0.5 / 0.0353,
                    legend = TRUE)

ggsave("./output_files/FigE5a_XO_dist_heat.pdf", XO_heat, dpi = 300,
       useDingbats=FALSE)     
