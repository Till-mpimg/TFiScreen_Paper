#!/bin/bash
# Master script to analyze identified TF-RE interactions and plot Fig. 4e, h-i, Extended Fig. E8c-j

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Plot total number of positive and negative interactions for Fig. 4e
echo -e "Rscript ${script_dir}Fig4e_plot_REPi_interactions.R $wd"
Rscript ${script_dir}Fig4e_plot_REPi_interactions.R $wd

#Plot numbers of interactions per TF/RE for Extended Fig. E8c-d and FACS data for Extended Fig. E8f
echo -e "Rscript ${script_dir}FigE8cdf_plot_numbers.R $wd"
Rscript ${script_dir}FigE8cdf_plot_numbers.R $wd

#Calculate correlation between reporter screens and TFi screen for significant interactions and plot Extended Fig. E8g
echo -e "Rscript ${script_dir}FigE8g_REPi_TFi_cor.R $wd"
Rscript ${script_dir}FigE8g_REPi_TFi_cor.R $wd

#Compare REPi to TFi screen data and plot Fig. 4h
echo -e "Rscript ${script_dir}Fig4h_REPI_cumulative_interaction.R $wd"
Rscript ${script_dir}Fig4h_REPI_cumulative_interaction.R $wd

#Plot expression for reporter screen interactors in Extended Fig. E8h-j
echo -e "Rscript ${script_dir}FigE8hij_REPi_expression.R $wd"
Rscript ${script_dir}FigE8hij_REPi_expression.R $wd

#Perform GSEA and plot Fig. 4i
echo -e "Rscript ${script_dir}Fig4i_REPi_gsea.R $wd"
Rscript ${script_dir}Fig4i_REPi_gsea.R $wd

#Perform motif analysis from vierstra data
echo -e "${script_dir}find_motifs.sh $wd"
${script_dir}find_motifs.sh $wd

#Clean FIMO results
echo -e "Rscript ${script_dir}clean_fimo.R $wd"
Rscript ${script_dir}clean_fimo.R $wd

#Quantify motif enrichment in strong interactions to plot Extended Fig. E8e
echo -e "Rscript ${script_dir}FigE8e_plot_fimo.R $wd"
Rscript ${script_dir}FigE8e_plot_fimo.R $wd
