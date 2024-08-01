library(tidyverse)
library(viridis)


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

#Read results
bulk_pca <- read.delim("./output_files/bulk_pca.tsv")
sc_pca <- read.delim("./output_files/pseudobulk_pca.tsv")

bulk_df <- bulk_pca %>% 
  select(PC1 = principal.component.0, PC2 = principal.component.1, sample = targets) %>% 
  separate(sample, c("sex", "rep", "timepoint"))

sc_df <- sc_pca %>% 
  select(PC1 = principal.component.0, PC2 = principal.component.1, sample = targets) %>% 
  mutate(sample = str_remove(sample, "clust_"))


plot <- bulk_df %>% 
  ggplot() +
  geom_point(aes(x = PC1, y = PC2, color = timepoint, shape = sex), size = 1.2) +
  geom_point(data = sc_df, aes(x = PC1, y = PC2, fill = sample), size = 2.5, shape = 21) +
  scale_color_viridis(option = "plasma", discrete = TRUE, end = 0.9) +
  scale_fill_viridis(option = "plasma", discrete = TRUE) +
  scale_shape_manual(values = c(2,1))

fix <- set_panel_size(plot, height = unit(4, "cm"), width = unit(4, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig3b_pca_bulk_sc.pdf", fix, dpi = 300,
       useDingbats=FALSE)
