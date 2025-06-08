#Script to plot manually counted RNA-FISH data
library(tidyverse)
library(egg)
library(gridExtra)


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

#Reads manual count tables (CZI files provided at doi.org/10.5281/zenodo.12821363)
fish_total <- read.delim("./input_files/FISH_sgOct4_manual_count.txt")
fish_sum <- fish_total %>% 
  mutate(xist_ex <- factor(xist_ex), rep = paste0("R", rep))  %>% 
  group_by(Med, sgRNA, rep, xist_ex) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = xist_ex, values_from = n) %>% 
  mutate(`1` = replace(`1`, is.na(`1`), 0), `2` = replace(`2`, is.na(`2`), 0), xist_perc = `1` + `2`)
  

#Perform t-test for %Xist+ cells between knockdown and control
test_df <- data.frame(Med = unique(fish_sum$Med))

test_df$pval <- mapply(function(a) {t.test(fish_sum[fish_sum$sgRNA=="sgNT" & fish_sum$Med==a,]$xist_perc,
                                             fish_sum[fish_sum$sgRNA=="sgOct4" & fish_sum$Med==a,]$xist_perc, 
                                             var.equal = TRUE)$p.value},
                       test_df$Med)

#Plot percentage of Xist positive cells
xist_plot <- fish_sum %>%
  ggplot(aes(x = factor(Med, levels = c("2i", "-2iL", "EPI")), y = xist_perc)) +
  stat_summary(aes(fill = factor(sgRNA, levels = c("sgNT", "sgOct4"))), 
                   geom = "bar", fun = "mean", position = position_dodge(width = 1), width = 0.8) +
  stat_summary(aes(group = factor(sgRNA, levels = c("sgNT", "sgOct4")), y = `2`), 
               geom = "bar", fun = "mean", position = position_dodge(width = 1), color = "black", fill = NA, width = 0.8) +
  geom_point(aes(fill = factor(sgRNA, levels = c("sgNT", "sgOct4"))),
             shape = 21, size = 0.7, position = position_dodge(width = 1)) +
  scale_y_continuous(limits = c(0, 100), name = "%Xist-positive cells") +
  scale_x_discrete(name = "Medium") +
  geom_text(data = na.omit(test_df[test_df$pval<=0.05,]), aes(label = "*", y = 85), 
            size = 8/2.8) +
  scale_fill_manual(values = c("#676767", "#F7A71C"), name = "Guide")

fix <- set_panel_size(xist_plot, height = unit(2, "cm"), width = unit(1.5, "cm"))
grid.arrange(fix)
ggsave("./output_files/Fig2c_Oct4_FISH_perc.pdf", fix, dpi = 300,
       useDingbats=FALSE)


