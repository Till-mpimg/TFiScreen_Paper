#!/bin/bash
# Master script to generate count tables and perform target enrichment for the TFi screen. Also performs quality control and plots Extended Fig. E1e-j

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Performs counting and target enrichment using MAGeCK
echo -e "${script_dir}TFi_mageck.sh $wd"
${script_dir}TFi_mageck.sh $wd

#Processes FACS data and plots Extended Fig. E1f
echo -e "Rscript ${script_dir}FigE1f_TFi_FACS.R $wd"
Rscript ${script_dir}FigE1f_TFi_FACS.R $wd

#Plot raw figures for Extended Fig. E1g
echo -e "Rscript ${script_dir}FigE1g_TFi_sorting_strategy.R $wd"
Rscript ${script_dir}FigE1g_TFi_sorting_strategy.R $wd

#Plot RNA-FISH data for Extended Fig. E2a
echo -e "Rscript ${script_dir}FigE2a_TFi_FISH.R $wd"
Rscript ${script_dir}FigE2a_TFi_FISH.R $wd

#Plot qPCR data for Extended Fig. E2b
echo -e "Rscript ${script_dir}FigE2b_TFi_qPCR.R $wd"
Rscript ${script_dir}FigE2b_TFi_qPCR.R $wd

#Visualize QC for Extended Fig. E2c-d
echo -e "Rscript ${script_dir}FigE2c-d_TFi_QC.R $wd"
Rscript ${script_dir}FigE2c-d_TFi_QC.R $wd

#Visualize target enrichment between Plasmid library and the Selected population. Plots Extended Fig. E2e
echo -e "Rscript ${script_dir}FigE2e_TFi_Selected_rank.R $wd"
Rscript ${script_dir}FigE2e_TFi_Selected_rank.R $wd
