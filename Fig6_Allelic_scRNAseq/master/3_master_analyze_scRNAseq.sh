#!/bin/bash
# Master script to analyze scRNAseq for chrX allelic expression

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Filters for monoallelic Xist-positive and plots Fig. 6c-d and Extended Fig. E10h
echo -e "Rscript ${script_dir}Fig6cdE10h_generate_filtered_table.R $wd"
Rscript ${script_dir}Fig6cdE10h_generate_filtered_table.R $wd

#Analyzes XCI in groups of cells with variable Xist levels. Creates gene clusters of the DEL scRNAseq data.
echo -e "Rscript ${script_dir}Cluster_X_genes_DEL.R $wd"
Rscript ${script_dir}Cluster_X_genes_DEL.R $wd

#Compares XCI between deletion and wildtype cells. Plots Fig. 6e and Extended Fig. E10i
echo -e "Rscript ${script_dir}Fig6eE10i_comp_DEL.R $wd"
Rscript ${script_dir}Fig6eE10i_comp_DEL.R $wd



