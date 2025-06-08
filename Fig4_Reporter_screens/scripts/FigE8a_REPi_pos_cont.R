# Analyzes fold change of positive controls in reporter screens and plots Extended Fig. E9a
library(tidyverse)
library(pheatmap)
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

counts <- read.delim("./input_files/REPi.count_normalized.txt")
lfc <- read.delim("./input_files/REPi_lfc.txt")
level <- c("RE57L", "RE57M", "RE57R", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127")

#Splits up RE57/58 guides according to location
sgRNA <- counts[counts$Gene == "RE57",]$sgRNA
RE58 <- c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
RE57R <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE)
RE57M <- c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)
RE57L <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
Rep_Gene <- c("RE58_RE57R", "RE57R", "RE57R_RE57M", "RE58", "RE57M", "RE57M", "RE57L", "RE58_RE57R", "RE57R", "RE57M")

RE57_df <- data.frame(sgRNA, Rep_Gene)

##Plot guide controls
cont_df <- lfc %>% 
  filter(Gene %in% c("RE57", "RE58", "RE61", "RE85", "RE93", "RE95", "RE96", "RE97", "RE127", "FIREWACh")) %>% 
  left_join(RE57_df) %>% 
  pivot_longer(c(Gene, Rep_Gene), names_to = "trash", values_to = "Gene") %>% 
  na.omit() %>% 
  filter(!Gene == "RE57") %>% 
  separate(Gene, c("Gene1", "Gene2")) %>% 
  pivot_longer(c(Gene1, Gene2), names_to = "trash2", values_to = "Gene") %>% 
  na.omit() %>% 
  filter(!(sgRNA == "gRNA_2528" & reporter != Gene & reporter %in% c("RE58", "RE57R"))) %>% 
  filter(!(sgRNA == "gRNA_2455" & reporter != Gene & reporter %in% c("RE57M", "RE57R"))) %>% 
  filter(!(sgRNA == "gRNA_2537" & reporter != Gene & reporter %in% c("RE58", "RE57R"))) 

cont_pheat <- cont_df %>% 
  group_by(Gene, reporter) %>% 
  summarize(lfc = mean(lfc)) %>% 
  select(Gene, reporter, lfc) %>% 
  pivot_wider(names_from = Gene, values_from = lfc) %>% 
  column_to_rownames("reporter") %>% 
  as.matrix()

mat_col <- colorRampPalette(c("#666666", "#FFFFFF"))(6)

cont_heat <- pheatmap(cont_pheat[c("noRE", level),c("FIREWACh", level)], 
                      border_color = "#C6C6C6",
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      fontsize = 6, 
                      color = mat_col, 
                      breaks = seq(-6, 0, 1),
                      cellwidth = 0.25 / 0.0353, 
                      cellheight = 0.25 / 0.0353,
                      legend = TRUE)

ggsave("./output_files/FigE8a_REPi_cont_heatmap.pdf",cont_heat, dpi = 300,
       useDingbats=FALSE)