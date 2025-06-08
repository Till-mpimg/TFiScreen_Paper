#!/bin/bash
# Master script to analyze RNA-FISH/qPCR from differentiation experiments and Oct4 knockdowns. Plots Fig. 2b-c and Extended Fig. E5a-d

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Analyzes and plots qPCR data for Fig. 2b and Extended Fig. E4c-d
echo -e "Rscript ${script_dir}Fig2bE4c-d_sgOct4_qPCR.R $wd"
Rscript ${script_dir}Fig2bE4c-d_sgOct4_qPCR.R $wd

#Plots manually counted RNA-FISH data following Oct4 knockdown as Fig. 2c
echo -e "Rscript ${script_dir}Fig2c_plot_FISH_manual_count.R $wd"
Rscript ${script_dir}Fig2c_plot_FISH_manual_count.R $wd

#Plots qPCR data for Extended Fig. E4a
echo -e "Rscript ${script_dir}FigE4a_-2iL_qPCR.R $wd"
Rscript ${script_dir}FigE4a_-2iL_qPCR.R $wd

#Analyzes and plots qPCR data for Extended Fig. E4b
echo -e "Rscript ${script_dir}FigE4b_EpiLC_qPCR.R $wd"
Rscript ${script_dir}FigE4b_EpiLC_qPCR.R $wd



