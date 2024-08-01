#!/bin/bash
# Master script to analyze RNA-seq time course and plot Fig. 3e-i and Extended Fig. E6b-e.

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Plots z-score expression dynamics for different gene groups and calculates Fisher enrichment to plot Fig. 3d-e, g-h and Extended Fig. E6b
echo -e "Rscript ${script_dir}Fig3deghE6b_plot_zscores.R $wd"
Rscript ${script_dir}Fig3deghE6b_plot_zscores.R $wd

#Plots expression dynamics of individual genes for Fig. 3c,i and Extended Fig. E6d
echo -e "Rscript ${script_dir}Fig3ciE6d_plot_indiv_genes.R $wd"
Rscript ${script_dir}Fig3ciE6d_plot_indiv_genes.R $wd

#Calculates ANOVA and plots Fig. 3f/Extended Fig. E6c
echo -e "Rscript ${script_dir}Fig3fE6c_calc_ANOVA.R $wd"
Rscript ${script_dir}Fig3fE6c_calc_ANOVA.R $wd

#Analyzes Pacini data and plot Extended Fig. E6e
echo -e "Rscript ${script_dir}FigE6e_plot_pacini_genes.R $wd"
Rscript ${script_dir}FigE6e_plot_pacini_genes.R $wd



