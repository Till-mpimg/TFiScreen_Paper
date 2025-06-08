#qpcr for Differentiation markers following Xist activator knockdown
library(tidyverse)
library(scales)
library(egg)
library(readxl)
library(EnvStats)

theme_set(theme_classic() + 
            theme(legend.title = element_blank(), 
                  axis.title.x = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), 
                  axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1), axis.title.y = element_text(size = 6), 
                  legend.text = element_text(size = 6), axis.text.y = element_text(size = 6),
                  strip.background = element_blank(), 
                  strip.text = element_text(size = 6)))


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Takes the standard .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC and removes all samples that did not include all target genes (e.g. negative controls).
raw <- list.files(path = "./input_files/", pattern = ".*DIFF.*.xls", full.names = TRUE) %>% 
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

write.table(output, "/project/ag_schulz/Till/TFi_Paper/revision/output/qPCR_relExp_DIFFmarkers.txt",sep="\t",row.names = FALSE, quote = FALSE)

plot_df <- output %>%  
  separate(Sample, c("rep", "med", "guide"), sep = "_") %>% 
  pivot_longer(-c(1:3), names_to = "Gene", values_to = "Rel_Exp") %>% 
  na.omit() 

#Calc significance
test_df <- expand.grid(Gene = unique(plot_df$Gene), guide = unique(plot_df$guide))

test_df$pval <- mapply(function(a,b)
{t.test(plot_df[plot_df$Gene==a & plot_df$guide=="LR15",]$Rel_Exp,
        plot_df[plot_df$Gene==a & plot_df$guide==b,]$Rel_Exp, 
        var.equal = TRUE)$p.value},
test_df$Gene, test_df$guide)

help_df <- plot_df %>% 
  group_by(Gene) %>% 
  summarize(pos = max(Rel_Exp + 0.5)) %>% 
  left_join(test_df)

#Plot DIFF markers
xist_plot <- plot_df %>% 
  ggplot(aes(x = factor(guide, levels = c("LR15", "TS18", "TS20", "TS16", "TS22")), 
             y = Rel_Exp, color = factor(guide, levels = c("LR15", "TS18", "TS20", "TS16", "TS22")))) +
  facet_wrap(~Gene, scales = "free_y") +
  geom_point(position = position_dodge(width = 1)) + 
  stat_summary(geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, color = "black", position = position_dodge(width = 1)) +
  ylab("Relative expression (log2)") + 
  geom_text(data = help_df[help_df$pval<=0.05,], aes(label = "*", y = pos), color = "black", 
            size = 8/2.8, position = position_dodge(width = 0.5)) + 
  scale_y_continuous(expand = c(0,1)) 

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE9h_DIFFmarkers_qpcr.pdf", fix, dpi = 300,
       useDingbats=FALSE)



