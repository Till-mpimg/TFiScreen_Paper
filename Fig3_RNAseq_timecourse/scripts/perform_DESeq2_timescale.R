library(DESeq2)
library(tidyverse)


#Specify working directory
args <- commandArgs(trailingOnly = TRUE)
wd <-  args[1]

setwd(wd)

#Reads RNA-seq timecourse counts
bulk_rnaseq_counts = read.table("./input_files/Counts_RNA_timecourse.txt", sep="\t", header = T)
rownames(bulk_rnaseq_counts)= bulk_rnaseq_counts[,2]

#Creates dataframe to map gene to gene_id
gene_sym_ens_id_mapping = bulk_rnaseq_counts[c("gene","gene_id")]

#Removes gene columns
bulk_rnaseq_counts = bulk_rnaseq_counts[,-c(1,2)]

#Prepares coldata for DESeq2
coldata = data.frame(colnames(bulk_rnaseq_counts))
colnames(coldata) = c("samples")
coldata = coldata %>% 
   separate(samples, into = c("sex", "rep", "timepoint"), sep="_",remove =F)
rownames(coldata) = coldata[["samples"]]
coldata$timepoint <- factor(coldata$timepoint)
coldata$sex <- factor(coldata$sex)
coldata$rep <- factor(coldata$rep)

bulk_rnaseq_counts <- bulk_rnaseq_counts[, rownames(coldata)]

#Runs DESeq2 for all timepoints
dds <- DESeqDataSetFromMatrix(countData = bulk_rnaseq_counts,
                              colData = coldata,
                              design = (~ timepoint  ))

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

#Return significant genes per timepoint (FDR<0.05, log2FC>1.5)
res_timepoint_sig_genes = c()
for (comparison in resultsNames(dds)[2:length(resultsNames(dds))]){
    res_timepoint <- results(dds, name=comparison, test="Wald",alpha = 0.05)
    res_timepoint_sig <- subset(res_timepoint, padj < 0.05)
    res_timepoint_sig <- res_timepoint_sig[ which( res_timepoint_sig$log2FoldChange > 1.5 | res_timepoint_sig$log2FoldChange < -1.5) , ]
    
    
    res_timepoint_sig_genes = c(res_timepoint_sig_genes, rownames(res_timepoint_sig))
}


#get the gene symbols for these genes
write.csv(gene_sym_ens_id_mapping[(unique(res_timepoint_sig_genes)),]$gene, "./data/alldiffGenes_timeseries.csv")

