#!/bin/bash
# Master script to analyze RNA-FISH/qPCR from differentiation experiments and Oct4 knockdowns. Plots Fig. 2b-c and Extended Fig. E5a-d

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Analyzes and plots qPCR data for Fig. 2b and Extended Fig. E5c-d
echo -e "Rscript ${script_dir}Fig2bE5c-d_sgOct4_qPCR.R $wd"
Rscript ${script_dir}Fig2bE5c-d_sgOct4_qPCR.R $wd

#Plots manually counted RNA-FISH data following Oct4 knockdown as Fig. 2b
echo -e "Rscript ${script_dir}Fig2c_plot_FISH_manual_count.R $wd"
Rscript ${script_dir}Fig2c_plot_FISH_manual_count.R $wd

#Plots qPCR data for Extended Fig. E5a
echo -e "Rscript ${script_dir}FigE5a_RA_qPCR.R $wd"
Rscript ${script_dir}FigE5a_RA_qPCR.R $wd

#Analyzes and plots qPCR data for Extended Fig. E5b
echo -e "Rscript ${script_dir}FigE5b_EpiLC_qPCR.R $wd"
Rscript ${script_dir}FigE5b_EpiLC_qPCR.R $wd



