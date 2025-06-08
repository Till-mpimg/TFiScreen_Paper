#!/bin/bash
# Master script to analyze RNA-seq time course and plot Fig. 3e-i and Extended Fig. E5b-f.

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'


#Plots z-score expression dynamics for different gene groups and calculates Fisher enrichment to plot Fig. 3d-e, g-h and Extended Fig. E5b
echo -e "Rscript ${script_dir}Fig3deghE5b_plot_zscores.R $wd"
Rscript ${script_dir}Fig3deghE5b_plot_zscores.R $wd

#Plots expression dynamics of individual genes for Fig. 3c,i and Extended Fig. E5d
echo -e "Rscript ${script_dir}Fig3ciE5d_plot_indiv_genes.R $wd"
Rscript ${script_dir}Fig3ciE5d_plot_indiv_genes.R $wd

#Calculates ANOVA and plots Fig. 3f/Extended Fig. E5c
echo -e "Rscript ${script_dir}Fig3fE5c_calc_ANOVA.R $wd"
Rscript ${script_dir}Fig3fE5c_calc_ANOVA.R $wd

#Analyzes Pacini data and plot Extended Fig. E5e
echo -e "Rscript ${script_dir}FigE5e_plot_pacini_genes.R $wd"
Rscript ${script_dir}FigE5e_plot_pacini_genes.R $wd

#Analyzes pseudobulk scRNAseq in vivo data and plot Extended Fig. E5f
echo -e "Rscript ${script_dir}FigE5f_Mohammed_zscore.R $wd"
Rscript ${script_dir}FigE5f_Mohammed_zscore.R $wd
