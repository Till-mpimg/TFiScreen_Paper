#qpcr for sgFtx-Xert-REs control line
library(tidyverse)
library(scales)
library(egg)
library(readxl)
library(EnvStats)

theme_set(theme_classic() + 
            theme(legend.title = element_blank(), 
                  axis.title.x = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), 
                  axis.text.x = element_text(size = 6), axis.title.y = element_text(size = 6), 
                  legend.text = element_text(size = 6), axis.text.y = element_text(size = 6),
                  strip.background = element_blank(), 
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Takes the standard .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC and removes all samples that did not include all target genes (e.g. negative controls).
raw <- read_xls("./input_files/qPCR_CRISPRi_Cont_Plate1.xls", sheet = "Results", skip = 40) %>% 
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

write.table(output, "./output_files/qPCR_CRISPRi_cont.txt",sep="\t",row.names = FALSE, quote = FALSE)

plot_df <- output %>%  
  separate(Sample, c("rep", "day", "line"), sep = "_") %>% 
  pivot_longer(-c(1:3), names_to = "Gene", values_to = "Rel_Exp") %>% 
  na.omit() 

#Calc significance
test_df <- expand.grid(Gene = unique(plot_df$Gene), cont = c("A3", "LR15")) %>%  
  mutate(line = ifelse(cont == "A3", "F9", "AMC7"), comp = ifelse(cont == "A3", "DEL", "SP107"))

test_df$pval <- mapply(function(a,b,c)
{t.test(plot_df[plot_df$Gene==a & plot_df$line==b,]$Rel_Exp,
        plot_df[plot_df$Gene==a & plot_df$line==c,]$Rel_Exp, 
        var.equal = TRUE)$p.value},
test_df$Gene, test_df$cont, test_df$line)


#Plot genes
xist_plot <- plot_df %>% 
  ggplot(aes(x = factor(line, levels= c("A3", "F9", "LR15", "AMC7")),
             y = Rel_Exp, color = factor(line, levels= c("A3", "F9", "LR15", "AMC7")))) +
  facet_wrap(~Gene, scales = "free_y") +
  geom_point(position = position_dodge(width = 1)) + 
  stat_summary(aes(group = line), geom = "crossbar", fun = "mean", 
               width = 0.5, lwd = 0.25, color = "black", position = position_dodge(width = 1)) +
  ylab("Relative expression (log2)") + 
  scale_color_manual(values = c("#121212", "#AA0023", "#121212", "#AA0023")) +
  scale_y_continuous(expand = c(0.2,0)) 

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("/project/ag_schulz/Till/TFi_Paper/revision/output/F9A3_Cont_qPCR_relExp.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Calculate LFC
lfc_df <-  plot_df %>% 
  pivot_wider(names_from = line, values_from = Rel_Exp) %>% 
  mutate(SP107 = AMC7 - LR15, DEL = F9 - A3) %>% 
  pivot_longer(c(SP107, DEL), values_to = "lfc", names_to = "comp")

dummy_data <- lfc_df %>%
  group_by(Gene) %>%
  summarize(comp = "DEL", lfc = c(-3, 3))

lfc_plot <- lfc_df %>% 
  ggplot(aes(x= comp, y = lfc)) +
  facet_wrap(~Gene, scales = "free") +
  geom_point(size = 0.8, color = "#676767") +
  geom_point(data = dummy_data, aes(x = comp, y = lfc), color = "white", alpha = 0) +  # Invisible points
  stat_summary(geom = "crossbar", fun = "mean", lwd = 0.25, width = 0.5)+
  ylab("LFC (WT/Perturbation)")  + 
  geom_text(data = test_df[test_df$pval<=0.05,], aes(label = "*", y = 2.5), color = "black", 
            size = 8/2.8, position = position_dodge(width = 0.5)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed")

fix <- set_panel_size(lfc_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE10b_qPCR_CRISPRi_cont.pdf", fix, dpi = 300,
       useDingbats=FALSE)
