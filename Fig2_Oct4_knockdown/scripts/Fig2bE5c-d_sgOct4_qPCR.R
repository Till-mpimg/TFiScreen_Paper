#qpcr
library(tidyverse)
library(scales)
library(egg)
library(readxl)
library(EnvStats)
library(gridExtra)
library(ggh4x)

theme_set(theme_classic() +
            theme(legend.title = element_blank(), 
                  axis.title.x = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), 
                  axis.text = element_text(size = 6), axis.title.y = element_text(size = 6), 
                  legend.text = element_text(size = 6), 
                  strip.background = element_blank(), 
                  strip.text = element_text(size = 6)))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)


#Takes the standard .xls file from the QPCR machine and puts the data into a wide format. Also removes all wells termed
#NTC and removes RT- samples.
raw <- list.files(path = "./input_files/", pattern = ".*sgOct4.*.xls", full.names = TRUE) %>% 
  lapply(read_xls, sheet = "Results", skip = 40) %>% 
  lapply(mutate, CT = ifelse(CT == "Undetermined", 40, CT)) %>%  
  lapply(transform, CT = as.numeric(CT)) %>% 
  bind_rows %>%
  select(Sample.Name, Target.Name, CT) %>%
  dplyr::rename(Sample = Sample.Name, Gene = Target.Name) %>%
  filter(Sample != "NTC") %>%
  filter(Gene != "Arpo-") %>% 
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarize(CT = geoMean(CT)) %>%
  na.omit()

analysis <- raw %>%
  pivot_wider(names_from = Gene, values_from = CT) %>% 
  mutate(Cont = (Arpo + Rrm2) /2, Arpo = NULL, Rrm2 = NULL) %>%
  select(Sample, Cont, everything()) 

#Performs ddCT analysis vs Arpo/Rrm2
for( i in 3:length(analysis) ) {
  analysis[i] <- analysis[i] - analysis[2]
  analysis[i] <- 2^ - analysis[i]
  analysis[i] <- log2(analysis[i])
}

output <- analysis %>%
  select(-Cont) 

#Writes relative expressions to .txt file
write.table(output, "./data/sgOct4_qPCR.txt",sep="\t",row.names = FALSE, quote = FALSE)

plot_df <- output %>% 
  separate(Sample, c("rep", "day", "med", "dTAG", "guide"), sep = "_") %>% 
  pivot_longer(-c(1:5), names_to = "Gene", values_to = "Rel_Exp") %>% 
  na.omit()


#Performs t-test comparing control and knockdown conditions
test_df_d0 <- expand.grid(guide = unique(plot_df$guide), Gene = unique(plot_df$Gene), day = "D0",
                          med = c("+", "2i"))

test_df_diff <- expand.grid(guide = unique(plot_df$guide), Gene = unique(plot_df[plot_df$day=="D2",]$Gene), day = "D2",
                            med = c("Epi", "DIFF", "RA"))

test_df <- rbind(test_df_d0, test_df_diff)

test_df$pval <- mapply(function(a,b,c)
{t.test(plot_df[plot_df$guide==a & plot_df$Gene==b & plot_df$med==c & plot_df$dTAG=="dT",]$Rel_Exp,
        plot_df[plot_df$guide==a & plot_df$Gene==b & plot_df$med==c & plot_df$dTAG=="KD",]$Rel_Exp, 
        var.equal = TRUE)$p.value},
test_df$guide, test_df$Gene, test_df$med)

test_df$pval <- format(test_df$pval, scientific = FALSE)

#Calculating log fold changes
lfc_df <- plot_df %>% 
  pivot_wider(names_from = dTAG, values_from = Rel_Exp) %>% 
  mutate(lfc = KD - dT) 

#Plot d0 and d2 seperately
pos_d0_lfc_df <- lfc_df %>% 
  filter(day == "D0") %>% 
  group_by(day, med, Gene) %>%
  summarize(max = max(lfc), range = max(lfc) - min(lfc)) %>% 
  mutate(star_pos = max + 0.2 * range)

help_d0_lfc_df <- test_df %>% 
  filter(day == "D0") %>% 
  mutate(sig = ifelse(pval <= 0.05, "sig", "ns")) %>% 
  left_join(pos_d0_lfc_df)


y_list_d0 <- list(scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-3, 3)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-3, 3)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
               scale_y_continuous(expand = c(0,1), breaks = c(-5, 0, 5), limits = c(-7, 7)),
               scale_y_continuous(expand = c(0,1), breaks = c(-5, 0, 5), limits = c(-7, 7)))

lfc_d0_plot <- lfc_df %>% 
  filter(day == "D0") %>% 
  ggplot(aes(x = factor(med, levels = c("2i", "+")), y = lfc, fill = guide, group = guide)) +
  facet_wrap(~Gene, scales = "free_y") +
  geom_point(position = position_dodge(width = 0.5), shape = 21) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  stat_summary(geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.5),
               color = "black") +
  facetted_pos_scales(y = y_list_d0) +
  ylab("Fold Change KD/dTAG [log2]") + 
  scale_fill_manual(values = c("#666666", "#F7A81F"))+
  geom_text(data = help_d0_lfc_df, aes(label = "*", y = star_pos, group = guide, color = sig), 
            size = 8/2.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c( "white", "black"))

fix <- set_panel_size(lfc_d0_plot, height = unit(2, "cm"), width = unit(1, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig2bE5c_sgOct4_qPCR_d0.pdf", fix, dpi = 300,
       useDingbats=FALSE)


#Plot day 2
pos_diff_lfc_df <- lfc_df %>% 
  filter(day != "D0") %>% 
  group_by(day, med, Gene) %>%
  summarize(max = max(lfc), range = max(lfc) - min(lfc)) %>% 
  mutate(star_pos = max + 0.2 * range)

help_diff_lfc_df <- test_df %>% 
  filter(day != "D0") %>% 
  mutate(sig = ifelse(pval <= 0.05, "sig", "ns")) %>% 
  left_join(pos_diff_lfc_df)


y_list_diff <- list(scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-3, 0, 3), limits = c(-5, 5)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-3, 3)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2, 2)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2.5, 2.5)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-2, 0, 2), limits = c(-2.5, 2.5)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-3, 0, 3), limits = c(-5, 5)),
                  scale_y_continuous(expand = c(0,1), breaks = c(-3, 0, 3), limits = c(-5, 5)))

lfc_diff_plot <- lfc_df %>% 
  filter(day != "D0") %>% 
  ggplot(aes(x = factor(med, levels = c("DIFF", "Epi", "RA")), y = lfc, fill = guide, group = guide)) +
  facet_wrap(day~Gene, scales = "free_y") +
  geom_point(position = position_dodge(width = 0.5), shape = 21) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  stat_summary(geom = "crossbar", fun = "mean", width = 0.5, lwd = 0.25, position = position_dodge(width = 0.5),
               color = "black") +
  facetted_pos_scales(y = c(y_list_diff, y_list_diff)) +
  ylab("Fold Change KD/dTAG [log2]") + 
  scale_fill_manual(values = c("#666666", "#F7A81F"))+
  geom_text(data = help_diff_lfc_df, aes(label = "*", y = star_pos, group = guide, color = sig), 
            size = 8/2.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c( "white", "black"))

fix <- set_panel_size(lfc_diff_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig2bE5d_sgOct4_qPCR_diff.pdf", fix, dpi = 300,
       useDingbats=FALSE)
