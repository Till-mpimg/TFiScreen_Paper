library(tidyverse)
library(egg)
library(gridExtra)
library(pheatmap)
library(viridis)
library(drc)
library(ggpubr)
library(pheatmap)

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]
setwd(wd)

xist_bin_col <- colorRampPalette(c("#ABABAB", "#93134D"))(10)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6), 
                  legend.title = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank()))

cell_df <- read.delim("./output_files/XistTitration_filtered_cell_table.txt")

#Bin cells according to Xist expression
bins <- cell_df %>% 
  dplyr::select(cell, Xist) %>% 
  unique() %>% 
  mutate(xist_bin = factor(ntile(Xist, n = 10)))

bins %>% group_by(xist_bin) %>% tally()

bin_df <- cell_df %>% 
  left_join(bins)
  
#Check if alleles are different
#Plot Xist expression and allelic ratio per bin
sum_bin_df <- bin_df %>%
filter(gene != "Xist") %>% 
group_by(Xist, cell, xist_bin) %>% 
summarize(sum_xi = sum(count_xi), sum_xa = sum(count_xa)) %>% 
mutate(chrX_ratio = sum_xi / (sum_xi + sum_xa))




xci_bin_plot <- sum_bin_df %>% 
  ggplot(aes(x = xist_bin, y = chrX_ratio, color = xist_bin)) +
  geom_violin(fill = NA) +
  geom_jitter(size = 0.2, alpha = 0.2, width = 0.2) +
  stat_summary(fun = "median", geom = "crossbar", color = "black", width = 0.5, lwd = 0.25) +
  scale_color_manual(values = xist_bin_col) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed")


allelic_fix <- set_panel_size(xci_bin_plot, height = unit(2, "cm"), width = unit(4, "cm"))
grid.arrange(allelic_fix)
ggsave("./output_files/Fig6k_XistTitration_chrX_ratio_bins.pdf", allelic_fix, 
       useDingbats=FALSE)



#Summarize allelic ratio/Bin
sum_genes <- bin_df %>% 
  group_by(xist_bin, gene, position) %>% 
  summarize(n = n(), ratio = mean(allelic_ratio))

#Filter genes with >=30 cells per bin
filter_genes <- sum_genes %>%
  ungroup() %>%  
  group_by(gene) %>% 
  summarize(min = min(n)) %>% 
  filter(min >= 30) %>% 
  dplyr::select(gene) %>% 
  unlist()

filter_sum <- sum_genes %>% 
  filter(gene %in% filter_genes)

out_table <- filter_sum %>% 
  dplyr::select(-position)

write_delim(out_table, "./output_files/XistTitration_bin_allelic_ratios.txt", delim = "\t")


#153 genes in the analysis
length(unique(filter_sum$gene))

#Normalize Xist expression to bin 10
xist_norm_df <- sum_bin_df %>% 
  mutate(Xist =  exp(Xist) - 1) %>% 
  group_by(xist_bin) %>% 
  summarize(Xist = median(Xist)) %>% 
  mutate(norm_Xist = Xist / Xist[10])

#CALC median/sd norm. xist in WT/DEL line
norm_factor = 18.6

old_data <- data.frame(sample = c("WT", "dFtx-Xert"), chrX_ratio = c(0.303, 0.4333), Xist = c(4.89, 2.15), 
                       Xist_low = c(2.33, 0.733), Xist_high = c(9.87, 4.32), 
                       low_ratio = c(0.183, 0.350), high_ratio = c(0.413, 0.477)) %>% 
  mutate(norm_Xist = Xist / norm_factor, norm_low = Xist_low / norm_factor, norm_high = Xist_high / norm_factor)

#Plot full chrX vs Xist
plot_bin_df <- sum_bin_df %>% 
  group_by(xist_bin) %>% 
  summarize(chrX_ratio = median(chrX_ratio)) %>% 
  left_join(xist_norm_df)

drm_model <- drm(chrX_ratio ~ norm_Xist, 
                 data = plot_bin_df, 
                 fct = LL.4(fixed = c(NA, 0, 0.5, NA)),
                 start = c(1, 0.25))

plot_bin_df$predicted <- predict(drm_model)

chrX_normXist <- plot_bin_df %>% 
  ggplot() +
  geom_line(aes(x = norm_Xist, y = predicted)) +
  geom_point(data = old_data, aes(x = norm_Xist, y = chrX_ratio), size = 2, color = "black") +
  geom_point(aes(x = norm_Xist, y = chrX_ratio), size = 1, color = "#93134D") +
  geom_segment(data = old_data, aes(x = norm_low, xend = norm_high, y = chrX_ratio, yend = chrX_ratio)) +
  geom_segment(data = old_data, aes(x = norm_Xist, xend = norm_Xist, y = low_ratio, yend = high_ratio)) +
  scale_x_continuous(limits = c(0, 1), name = "norm. Xist expression [to max bin]") +
  scale_y_continuous(limits = c(0, 0.55), name = "Allelic Ratio [chrX]", breaks = c(0, 0.25, 0.5))

model_fix <- set_panel_size(chrX_normXist, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(model_fix)
ggsave("./output_files/Fig6l_XistTitrration_full_chrX_doseresponse.pdf", model_fix, 
       useDingbats=FALSE)


#Calculate Hill Coefficients and ED50s
hill_fun <- function(x) {
  hill_df <- filter_sum %>% 
    filter(gene == x) %>% 
    left_join(xist_norm_df)
  
  min_target <- min(hill_df$ratio)
  max_target <- 0.5
  
  drm_model <- drm(ratio ~ norm_Xist, 
                   data = hill_df, 
                   fct = LL.4(fixed = c(NA, 0, max_target, NA)),
                   start = c(1, 0.25))
  
  #Get out model parameters
  ED50 <- ED(drm_model, 0.25, interval = "delta", type = "absolute")[1]
  hill_coef <- unname(drm_model$coefficients[1])
  
  #Calculate R squared to assess model fit
  tss <- sum((hill_df$ratio - mean(hill_df$ratio))^2)
  predicted_values <- predict(drm_model)
  rss <- sum((hill_df$ratio - predicted_values)^2)
  rsq <- 1 - (rss / tss)
  
  
  out_vec <- c(unique(hill_df$gene), hill_coef, ED50, min_target, rsq)
  names(out_vec) <- NULL
  
  return(list(model = drm_model, params = out_vec))
}

escape_genes <- filter_sum %>% 
  dplyr::select(-n) %>% 
  filter(xist_bin %in% c(1,10)) %>% 
  pivot_wider(names_from = xist_bin, values_from = ratio) %>% 
  mutate(total_silencing = `1` - `10`) %>% 
  filter(total_silencing <= 0.15)


run_vec <- unique(filter_sum$gene)
run_vec <- run_vec[!run_vec %in% escape_genes$gene]

out_list <- lapply(run_vec, function(x) hill_fun(x))

out_df <- as.data.frame(do.call(rbind, lapply(out_list, function(x) x$params))) %>% 
  filter(V2 <= 10)
colnames(out_df) <- c("gene", "hill_coef", "ED50", "min_ratio", "rsq")

out_text <- out_df %>% 
  dplyr::select(gene, hill_coef, ED50, min_ratio, rsq)

write_delim(out_text, "./output_files/XistTitration_dose_response_out.txt", delim = "\t")

#Plot rsquared as histogram
hist_plot <- out_df %>% 
  ggplot(aes(x = as.numeric(rsq))) +
  geom_histogram(binwidth = 0.05, fill = "#ABABAB") +
  scale_x_continuous(limits = c(0, 1), name = "R-squared") +
  geom_vline(aes(xintercept = 0.7), linetype = "dashed", color = "darkred")


hist_fix <- set_panel_size(hist_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(hist_fix)
ggsave("/project/ag_schulz/Till/TFi_Paper/revision/output/scXist_rsquared.pdf", hist_fix, 
       useDingbats=FALSE)

out_df <- out_df[out_df$rsq >= 0.7,]


#Plot the ED50s
ed50_median <- median(as.numeric(out_df$ED50))

ed50_plot <- out_df %>% 
  ggplot(aes(x = as.numeric(ED50))) +
  geom_histogram(binwidth = 0.05, fill = "#ABABAB") +
  scale_x_continuous(name = "ED50", limits = c(0, 1)) +
  geom_vline(aes(xintercept = ed50_median), color = "darkred", linetype = "dashed")


ed50_fix <- set_panel_size(ed50_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(ed50_fix)
ggsave("./output_files/Fig6n_XistTitration_doseresponse_ed50.pdf", ed50_fix, 
       useDingbats=FALSE)

#Relate model parameters to silencing speed / position / expression
#Silencing statistics (Barros de Sousa 2019)
sil_statistics <- read.delim("./input_files/sousa_silencing.txt") %>% 
  dplyr::select(gene = gene.name, half_time = half.time)

#Day 4 expression
xo_exp <- read.delim("./input_files/CPM_RNA_timecourse.txt") %>% 
  dplyr::select(gene, XO_R1_96h, XO_R2_96h, XO_R3_96h) %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "cpm") %>% 
  group_by(gene) %>% 
  summarize(xo_exp = mean(cpm))

#Add data to out_df
pos_df <- filter_sum %>% 
  ungroup() %>% 
  dplyr::select(gene, position) %>% 
  unique()

char_df <- out_df %>% 
  left_join(pos_df) %>% 
  left_join(sil_statistics) %>% 
  left_join(xo_exp) %>% 
  mutate(lin_dist = abs(position - 103483254)) %>%
  dplyr::select(-rsq) %>% 
  pivot_longer(c(hill_coef, ED50), names_to = "model_param", values_to = "value" )

#~out table
out <- char_df %>% 
  pivot_wider(names_from = model_param, values_from = value) 

write_delim(out, "./out_files/XistTitration_silencing_out.txt", delim = "\t")

#Plot Examples
plot_vec <- c("Rlim", "Zfp280c", "Wdr45", "Tfe3")


plot_fun <- function(gene_X) {
  gene_model <- Filter(function(x) x$params[1] == gene_X, out_list)[[1]]$model
  
  plot_df <- filter_sum %>% 
    filter(gene == gene_X) %>% 
    left_join(xist_norm_df)
  
  plot_df$predicted <- predict(gene_model)
  
  txt_df <- out_df %>% 
    filter(gene == gene_X)
  
  
  example_models <- plot_df %>% 
    ggplot(aes(x = norm_Xist)) +
    geom_line(aes(y = predicted)) +
    geom_point(aes(y = ratio), size = 1, color = "#93134D") +
    scale_x_continuous(limits = c(0, 1), name = "norm. Xist expression [to max bin]") +
    scale_y_continuous(limits = c(0, 0.55), name = "Allelic Ratio", breaks = c(0, 0.25, 0.5)) +
    geom_text(data = txt_df, aes(label = paste0("ED50: ", round(as.numeric(ED50), 2)), x = 0.7, y = 0.4), size = 6/2.8)+
    geom_text(data = txt_df, aes(label = paste0("n: ", round(as.numeric(hill_coef), 2)), x = 0.7, y = 0.5), size = 6/2.8)
  
  model_fix <- set_panel_size(example_models, height = unit(2, "cm"), width = unit(2, "cm"))
  grid.arrange(model_fix)
  ggsave(paste0("./output_files/Fig6o_XistTitration_doseresponse_", gene_X, ".pdf"), model_fix, 
         useDingbats=FALSE)
}

sapply(plot_vec, plot_fun)


#Investigate resitant genes on characteristics
kmeans_df <- out_df %>% 
  dplyr::select(gene, ED50) %>% 
  filter(ED50 <= 3) %>% 
  rownames_to_column("trash") %>% 
  dplyr::select(-trash) %>% 
  column_to_rownames("gene")


kmeans <- kmeans(kmeans_df, 3)
kmeans_df$cluster <- kmeans$cluster

kmeans_plot <- kmeans_df %>% 
  ggplot(aes(y = as.numeric(ED50), x = factor(cluster, levels = c("1", "3", "2")), fill = factor(cluster))) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_y_continuous(name = "ED50")

kmeans_fix <- set_panel_size(kmeans_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(kmeans_fix)
ggsave("./output_files/FigE10r_XistTitration_kmeans_genes.pdf", kmeans_fix, 
       useDingbats=FALSE)


type_df <- char_df %>% 
  left_join(kmeans_df) %>% 
  na.omit()


write_delim(type_df, "./output_files/XistTitration_kmeans.txt", delim = "\t")

exp_kmeans <- type_df %>% 
  ggplot(aes(x = factor(cluster, levels = c("1", "3", "2")), y = log10(xo_exp), fill = factor(cluster))) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(name = "Kmeans cluster")

kmeans_fix <- set_panel_size(exp_kmeans, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(kmeans_fix)
ggsave("./output_files/FigE10u_kmeans_exp_box.pdf", kmeans_fix, 
       useDingbats=FALSE)  

dist_kmeans <- type_df %>% 
  ggplot(aes(x = factor(cluster, levels = c("1", "3", "2")), y = lin_dist / 1000000, fill = factor(cluster))) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(name = "Kmeans cluster") +
  scale_y_continuous(name = 'Distance to Xist [Mb]')

kmeans_fix <- set_panel_size(dist_kmeans, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(kmeans_fix)
ggsave("./output_files/FigE10t_kmeans_dist_box.pdf", kmeans_fix, 
       useDingbats=FALSE)  

sil_kmeans <- type_df %>% 
  ggplot(aes(x = factor(cluster, levels = c("1", "3", "2")), y = half_time, fill = factor(cluster))) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(name = "Kmeans cluster") +
  scale_y_continuous(name = 'Silencing half-time [days]')

kmeans_fix <- set_panel_size(sil_kmeans, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(kmeans_fix)
ggsave("./output_files/FigE10s_kmeans_sil_box.pdf", kmeans_fix, 
       useDingbats=FALSE)  

#Plot total heatmap of the genes
order <- filter_sum %>% 
  group_by(gene) %>% 
  summarize(ratio = mean(ratio)) %>% 
  arrange(ratio)

pheat_mat <- filter_sum  %>%
  dplyr::select(-n, -position) %>%
  pivot_wider(names_from = xist_bin, values_from = ratio) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

dist_mat <- dist(pheat_mat)
hc <- hclust(dist_mat, method = "ward.D2")
plot(hc)

mat_col <- colorRampPalette(c("#ca0020", "#FFFFFF", "#0571b0"))(15)
mat_to_plot <- pheat_mat[order$gene, ]
rownames(mat_to_plot)[!(rownames(mat_to_plot) %in% c("Rlim", "Zfp280c", "Wdr45", "Tfe3"))] <- ""

target_heat <- pheatmap(mat_to_plot, 
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE,
                        border_color = NA,
                        show_rownames = TRUE,
                        fontsize = 6, 
                        color = mat_col, 
                        breaks = seq(0, 1.06666667, 0.06666667), 
                        cellwidth = 0.25 / 0.0353, 
                        cellheight = 0.05 / 0.0353, 
                        angle_col = 0)

ggsave("./output_files/Fig6m_XistTitration_heat.pdf", target_heat, 
       useDingbats=FALSE)  
