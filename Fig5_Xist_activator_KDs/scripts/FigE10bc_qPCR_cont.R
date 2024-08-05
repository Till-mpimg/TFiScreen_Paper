#qpcr cont for CUT&Tag/FISH experiments
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
raw <- list.files(path = "./input_files/", pattern = "Cont.*.xls", full.names = TRUE) %>% 
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

write.table(output, "./data/qPCR_Cont_relExp.txt",sep="\t",row.names = FALSE, quote = FALSE)

plot_df <- output %>%  
  separate(Sample, c("rep", "med", "guide"), sep = "_") %>% 
  pivot_longer(-c(1:3), names_to = "Gene", values_to = "Rel_Exp") %>% 
  na.omit()


#Plotting raw values
xist_plot <- plot_df %>% 
  filter(Gene == "Xist") %>%
  filter(med == "D2") %>% 
  ggplot(aes(x = factor(guide, levels = c("LR15", "TS18", "TS20", "TS16", "TS22")), 
             y = Rel_Exp, color = factor(guide, levels = c("LR15", "TS18", "TS20", "TS16", "TS22")))) +
  geom_point(position = position_dodge(width = 1)) + 
  stat_summary(geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, color = "black", position = position_dodge(width = 1)) +
  ylab("Relative expression (log2)") + 
  scale_y_continuous(expand = c(0,1)) 

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE10c_Cont_qPCR_Xist.pdf", fix, dpi = 300,
       useDingbats=FALSE)

#Calculate percent KD
help_kd <-  plot_df %>% 
  filter(Gene != "Xist") %>% 
  filter(guide == "LR15") %>% 
  select(rep, med, Gene, NT = Rel_Exp)

kd_lfc <- plot_df %>% 
  filter(Gene != "Xist") %>%  
  filter(guide != "LR15") %>% 
  left_join(help_kd) %>% 
  mutate(lfc = Rel_Exp - NT, perc = 2^Rel_Exp/(2^Rel_Exp + 2^NT)) 


kd_perc <- kd_lfc %>% 
  filter(med == "D2") %>% 
  ggplot(aes(x = factor(Gene, levels = c("Zic3", "Nfrkb", "Otx2", "Foxd3")),
             y = perc, color = factor(Gene, levels = c("Zic3", "Nfrkb", "Otx2", "Foxd3")))) +
  geom_point() + 
  stat_summary(geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, color = "black") +
  ylab("Percent expression of NT") + 
  scale_y_continuous(limits = c(0, 1)) 

fix <- set_panel_size(kd_perc, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE10b_Cont_qPCR_KDeff.pdf", fix, dpi = 300,
       useDingbats=FALSE)
