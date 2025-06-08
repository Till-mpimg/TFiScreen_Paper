library(tidyverse)
library(EnvStats)
library(Rgb)
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


transient1 <- c("Pou5f1", "Zic3", "Nfrkb", "Zfp518b", "Rbpj", "Sp1", "Sall4")
transient2 <- c("Otx2", "Nfe2l2", "Tox4", "Foxd3", "Zfp207", "Smarcc1", "Atf5", "Arid1a")
transient3 <- c("Epas1", "Myrf", "Mbd3", "Zscan10", "Usf1", "Arid1a", "Arnt")

basal <- c("Pou5f1", "Nfrkb", "Zic3", "Zfp518b", "Nfe2l2", "Adnp", "Tox4", "Zfp207", "Sp1", "Pias1", "Myrf", "Sall4", "Arid4a")
boost <- c("Otx2", "Rbpj", "Foxd3", "Arnt", "Zfp217", "Arid1a", "Hmg20b", "Epas1", "Mbd3", "Zscan10", "Usf1")

gencode_xert <- "./input_files/GENCODE_vM25_plus_Xert.gtf"
gene_names <- read.gtf(gencode_xert) %>%
  select(GENCODE = gene_id, gene = gene_name) %>%
  unique()

counts <- read.delim("./input_files/mohammed_2017_counts.txt") 

sex <- read.delim("./input_files/Mohammed_2017_sex_meta.txt")


cpm <- counts %>%
  mutate_at(vars(-GENCODE), funs( . / sum(.) * 10000))  %>%
  mutate_at(vars(-GENCODE), funs(log2(. + 1))) %>% 
  left_join(gene_names)


filt_cpm <- cpm %>% 
  filter(gene %in% c(transient1, transient2, transient3, basal, boost)) %>% 
  mutate(cluster = ifelse(gene %in% transient1, "transient1", 
                          ifelse(gene %in% transient2, "transient2",
                                 ifelse(gene %in% transient3, "transient3", "none")))) %>% 
  mutate(type = ifelse(gene %in% basal, "basal", 
                       ifelse(gene %in% boost, "boost", "none")))

filt_long <- filt_cpm %>% 
  pivot_longer(-c(GENCODE, gene, cluster, type), names_to = "cell", values_to = "cpm") %>% 
  left_join(sex) %>% 
  filter(stage != "E6.75")

filt_sum <- filt_long %>% 
  filter(cluster != "none") %>% 
  group_by(gene, stage, sex, cluster) %>% 
  summarize(mean_cpm = mean(cpm))

##Calculate based on zscore
z_sum <- filt_sum %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(zscore = scale(mean_cpm))

#t-test
test_df <- expand.grid(cluster = unique(z_sum$cluster), stage = unique(z_sum$stage))

test_fun <- function(a,b){
  t.test(z_sum[z_sum$cluster==a & z_sum$stage==b & z_sum$sex=="XX",]$zscore,
         z_sum[z_sum$cluster==a & z_sum$stage==b & z_sum$sex=="XY",]$zscore, paired = TRUE, 
         var.equal = TRUE)$p.value
}

test_df$pval <- mapply(test_fun, test_df$cluster, test_df$stage)
write_delim(test_df, "./data/FigE5f_Mohammed_zscore_pvals.txt", delim = "\t")

plot_groups <- z_sum  %>% 
  ggplot(aes(x = stage, y = zscore, color = sex)) +
  facet_wrap(~cluster) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = sex, group = sex), position = position_dodge(width = 0.8), size = 0.5) +
  scale_color_manual(values = c("#676767", "#CDCDCD"))


fix <- set_panel_size(plot_groups, height = unit(2, "cm"), width = unit(3, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE5f_Mohammed_zscore.pdf", fix, dpi = 300,
       useDingbats=FALSE)

