#This script was used to produce the TFi Library
library(tidyverse)
library(Rgb)
library(Rsubread)

theme_set(theme_classic() + 
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6),
                  strip.background = element_blank(), legend.title = element_blank()))

#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Specify source directory for TTseq data (doi.org/10.5281/zenodo.12822424, https://github.com/EddaSchulz/Xert_paper)
TTseq_bam <- "./input_files/ttseq_bam" 


#Adds transcription factors specified by AnimalTFdb 3.0 (Hu et al. 2019)
tfs <- read.delim("./input_files/Mus_musculus_TF_AnimalTFdb.txt") %>%
  select(gene_id = Ensembl)

#Adds the TPM counts from TX rnaseq timecourse (doi.org/10.5281/zenodo.12822424, Pacini et al. 2021)
RNA_TPM <- read.delim("./input_files/TPM_rnaseq_pacini.txt") %>%
  select(-contains("XO"))  %>% 
  mutate(gene_id = gsub("\\..*", "", gene_id))

#Adds GENCODE data supplemented with Xert coordinates (doi.org/10.5281/zenodo.12822424, Gjaltema et al. 2022)
transcripts_gencode <- read.gtf("./input_files/GENCODE_vM25_plus_Xert.gtf") %>%
  filter(feature == "transcript")

TSS_gencode <- transcripts_gencode %>% 
  mutate(TSS = ifelse(as.character(strand) == "+", start, end)) %>% 
  select(gene = gene_name, transcript = transcript_id, chromosome = seqname, strand, source, TSS)



#Adds putative X-linked activators from Ravid-Lustig et al. 2023 (FDR <= 0.05)
#Adds previously described Xist activators and binding proteins as controls
gene <- c("Xert", "Dcp1a", "Tsix", "Rif1", "Mbd1", "Brd4", "Spen", "Dnmt1", "Trim28", "Kat8", 
          "Eed", "Chd8", "Msl1", "Msl2", "Kansl1", "Kansl3", "Lbr", "Hnrnpu", "Hnrnpk")
Lustig_screen <- read.delim("./input_files/LustigScreen.gene_summary.txt") %>% 
  filter(Top.beta >= 0 & Top.wald.fdr <= 0.05) %>% 
  separate(Gene, c("id", "gene"), sep = "_")
controls <- data.frame(gene) %>%
  rbind(data.frame(gene = Lustig_screen$gene)) %>% 
  left_join(RNA_TPM)  %>% 
  transmute(gene_id = gene_id, gene = gene, XX_0d = (I_0h_XX + II_0h_XX + III_0h_XX) / 3,
            XX_1d = (I_24h_XX + II_24h_XX + III_24h_XX) / 3, XX_2d = (I_48h_XX + II_48h_XX + III_48h_XX) / 3,
            XX_3d = (I_72h_XX + II_72h_XX + III_72h_XX) / 3, XX_4d = (I_96h_XX + II_96h_XX + III_96h_XX) / 3)

#Adds the RNA-seq counts to the transcription factors and finds the maximum expression per factor
tf_counts <- tfs %>%
  inner_join(RNA_TPM) %>% 
  transmute(gene_id = gene_id, gene = gene, XX_0d = (I_0h_XX + II_0h_XX + III_0h_XX) / 3,
            XX_1d = (I_24h_XX + II_24h_XX + III_24h_XX) / 3, XX_2d = (I_48h_XX + II_48h_XX + III_48h_XX) / 3,
            XX_3d = (I_72h_XX + II_72h_XX + III_72h_XX) / 3, XX_4d = (I_96h_XX + II_96h_XX + III_96h_XX) / 3)
tf_counts$max <- do.call(pmax, tf_counts[,-c(1:2)])


#Plot max(TPM) distribution as a density plot. The Cut-off for expression used for the library design is indicated in yellow.
max_plot <- tf_counts %>%
  ggplot() +
  geom_density(aes(x = log2(max), ..scaled..)) +
  geom_vline(aes(xintercept = log2(10)), color = "#FAB336") +
  annotate(geom = "text", x = 0, y = 0.7, label = "Max TPM = 10", color = "#FAB336", size = 6/2.8) +
  annotate(geom = "text", x = 5, y = 0.2, label = "651/1636", color = "#FAB336", size = 6/2.8) +
  scale_x_continuous(name = "Max TPM (Log2)") +
  scale_y_continuous(name = "Guide Density [modal]", breaks = c(0, 0.5, 1))

fix <- set_panel_size(max_plot, height = unit(2, "cm"), width = unit(2, "cm"))
grid.arrange(fix)
ggsave("./output_files/FigE1b_maxTPM_distribution.pdf", fix, dpi = 300,
       useDingbats=FALSE)



#Filters the list on at least 10 TPM in any sample and adds the controls
tf_filtered <- tf_counts %>%
  filter(max >= 10) %>% 
  select(-max) %>% 
  rbind(controls) %>% 
  unique()

#Filters the transcript list on genes within the TF List and creates a sense and datasense dataframe for Featurecounts
#Uses TX TT-seq data to recover active isoforms (Gjaltema et al. 2022)
transcripts_gencode_filtered <- inner_join(TSS_gencode, tf_filtered, by = "gene") 
saf_up <- transcripts_gencode_filtered %>%
  transmute(GeneID = paste(gene, transcript, sep = "_"), Chr = chromosome, Start = ifelse(strand == "+", TSS - 2000, TSS),
            End = ifelse(strand == "+", TSS, TSS + 2000), Strand = strand)
saf_down <- transcripts_gencode_filtered %>%
  transmute(GeneID = paste(gene, transcript, sep = "_"), Chr = chromosome, Start = ifelse(strand == "+", TSS, TSS - 2000),
            End = ifelse(strand == "+", TSS + 2000, TSS), Strand = strand)

#FeatureCounts of 2000 bps up- and downstream of the annotated TSSs
#Files should be named by pattern: TT_TX1072_XX_d[0/2/4]_r[1/2]_dedup.bam
setwd(TTseq_bam)
temp = list.files(pattern="XX.*bam$")

feature_counts_up <- featureCounts(temp, annot.ext = saf_up, isPairedEnd = TRUE,
                                   strandSpecific = 2, allowMultiOverlap = TRUE)
counts_up <- data.frame(feature_counts_up$counts, feature_counts_up$annotation)
Dear Till Schwämmle,

feature_counts_down <- featureCounts(temp, annot.ext = saf_down, isPairedEnd = TRUE,
                                     strandSpecific = 2, allowMultiOverlap = TRUE)
counts_down <- data.frame(feature_counts_down$counts, feature_counts_down$annotation)


setwd(wd)


#Calculates log2 fold-change between the area upstream and downstream of the TSS and filters TSSs with LFC >= 1
up_down_lfc <- cbind(counts_down[7:12], log2((counts_down[1:6] + 0.01) / (counts_up[1:6] + 0.01))) %>% 
  mutate(TT_d0 = (.[[7]] + .[[8]]) / 2,
         TT_d2 = (.[[9]] + .[[10]]) / 2,
         TT_d4 = (.[[11]] + .[[12]]) / 2)
up_down_lfc$max <- do.call(pmax, up_down_lfc[,13:15]) 

filter_tss <- up_down_lfc %>% 
  filter(max >= 1) %>%
  separate(GeneID, c("gene", "transcript"), sep = "_")

transcripts_gencode_active <- transcripts_gencode_filtered %>%
  filter(transcript %in% filter_tss$transcript)


#Creates individual dataframes for all genes and puts them in a list
df_list <- transcripts_gencode_active %>%
  group_split(gene)


#Splits the df_list in 1 transcript and multiple transcripts. Assigns the single transcripts as cluster pA
df_list_mas <- df_list[sapply(df_list, function(x) nrow(x)[1]) > 1]
df_list_uno <- df_list[sapply(df_list, function(x) nrow(x)[1]) == 1]
df_list_single_transcript <- mapply(cbind, df_list_uno, cluster = 1, SIMPLIFY=F)

#Clusters the TSSs according to starting position with a height of 500 (should make bins ~500 bp big)
cluster_df_list <- lapply(df_list_mas, function(x) cutree(hclust(dist(x$TSS)), h = 500))

#Puts everything back in a single dataframe
df_list_clustered <- mapply(cbind, df_list_mas, cluster = cluster_df_list, SIMPLIFY=F)
df_list_final <- c(df_list_clustered, df_list_single_transcript)
df_final <- bind_rows(df_list_final) %>% 
  select(gene, transcript, chromosome, strand, source, TSS, cluster) %>%
  mutate(cluster = paste("p", cluster, sep = ""))  %>%
  mutate(id = paste(gene, cluster, sep = "_"))

write_delim(df_final, "./data/TFi_isoforms_raw.txt",
              delim = "\t", col_names = TRUE)

#Reduces the final dataframe to the TSS clusters(pA, pB, pC...). It puts start and end as the min/max TSS in the cluster.
df_reduced <- df_final %>% 
  group_by(gene, chromosome, strand, cluster, id) %>%
  summarize(start = min(TSS), end = max(TSS))


#Extends all promoters to 500
guidescan_df <- df_reduced %>% 
  mutate(length = end - start) %>% 
  mutate(start = ifelse(length <= 500 & strand == "-", start - (500 - length), start)) %>% 
  mutate(end = ifelse(length <= 500 & strand == "+", end + (500 - length), end)) %>% 
  mutate(length = end - start) %>% 
  mutate(Region.name = paste(chromosome, ":", start, "-", end, sep = ""))

write_delim(guidescan_df, "./TFi_target_promoters_500bp.txt", delim = "\t",
            col_names = TRUE)

#Creates bed file
guidescan_bed <- guidescan_df %>% 
  ungroup() %>% 
  select(chromosome, start, end, id) %>% 
  write_delim("./output_files/TFi_guidescan_500bp.bed", delim = "\t",
              col_names = FALSE)


##Guidescan_bed was transformed to mm39 via the UCSC LiftOver webtool 
##(at time of library design [2021/10/25] GuideScan2 only accepted mm39)
