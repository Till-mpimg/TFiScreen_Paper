#Script analyzing qPCRs from the validation of the TFi Screen
library(tidyverse)
library(scales)
library(egg)
library(readxl)
library(EnvStats)
library(gridExtra)
library(ggh4x)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), legend.text = element_text(size = 6),
                  axis.text = element_text(size = 6), axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
                  panel.border = element_rect(colour = "black", fill = NA, size = 0.5), axis.line = element_blank(), 
                  strip.background = element_blank(), axis.title.y = element_text(size = 6),
                  axis.title.x = element_blank(),
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Takes the standard .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC
raw <- list.files(path = "./input_files/", pattern = ".*EpiLC.*.xls", full.names = TRUE) %>% 
  lapply(read_xls, sheet = "Results", skip = 40) %>%  
  lapply(mutate, CT = ifelse(CT == "Undetermined", 40, CT)) %>%  
  lapply(transform, CT = as.numeric(CT)) %>% 
  bind_rows() %>%
  select(`Sample Name`, `Target Name`, "CT") %>%
  dplyr::rename(Sample = `Sample Name`, Gene = `Target Name`) %>%
  filter(Sample != "NTC") %>%
  filter(!str_detect(Sample, "-")) %>% 
  transform(CT = as.numeric(CT)) %>%
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarize(CT = geoMean(CT)) %>%
  na.omit()

#Calculate relative expression to Arpo and Rrm2
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

write.table(output, "./data/EpiLC_qPCR_RelExp.txt",sep="\t",row.names = FALSE, quote = FALSE)




qpcr <- output %>% 
  separate(Sample, c("rep", "day", "med"), sep = "_") %>% 
  pivot_longer(-c(rep, med, day), names_to = "gene", values_to = "rel_exp") %>% 
  na.omit()


#plot relative expression

y_list <- list(scale_y_continuous(expand = c(0,1), breaks = c(-9, -7, -5)),
               scale_y_continuous(expand = c(0,1), breaks = c(-5, -3, -1)),
               scale_y_continuous(expand = c(0,1), breaks = c(0, 1, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-5, -3, -1)),
               scale_y_continuous(expand = c(0,1), breaks = c(-16, -12, -8)),
               scale_y_continuous(expand = c(0,2), breaks = c(-10, -5, 0)))

plot <- qpcr %>% 
  ggplot(aes(x = med, y = rel_exp, color = med)) +
  facet_wrap(~gene, scales = "free") +
  geom_point(size = 0.7) +
  stat_summary(fun = "mean", geom = "crossbar", color = "black",
               lwd = 0.25) +
  scale_color_manual(values = c("#676767", "#F7A81F", "#F7A81F", "#F7A81F", "#F7A81F")) +
  facetted_pos_scales(y = y_list) 

fix <- set_panel_size(plot, height = unit(1.5, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5b_EpiLC_qPCR.pdf", fix, dpi = 300,
       useDingbats=FALSE)
