#!/bin/bash
# Master script to generate count tables and perform target enrichment for the TFi screen. Also performs quality control and plots Extended Fig. E1e-j

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Performs counting and target enrichment using MAGeCK
echo -e "${script_dir}TFi_mageck.sh $wd"
${script_dir}TFi_mageck.sh $wd

#Processes FACS data and plots Extended Fig. E1e
echo -e "Rscript ${script_dir}FigE1e_TFi_FACS.R $wd"
Rscript ${script_dir}FigE1e_TFi_FACS.R $wd

#Plot RNA-FISH data for Extended Fig. E1f
echo -e "Rscript ${script_dir}FigE1f_TFi_FISH.R $wd"
Rscript ${script_dir}FigE1f_TFi_FISH.R $wd

#Plot qPCR data for Extended Fig. E1g.
echo -e "Rscript ${script_dir}FigE1g_TFi_qPCR.R $wd"
Rscript ${script_dir}FigE1g_TFi_qPCR.R $wd

#Visualize QC for Extended Fig. E1h-i
echo -e "Rscript ${script_dir}FigE1h-i_TFi_QC.R $wd"
Rscript ${script_dir}FigE1h-i_TFi_QC.R $wd

#Visualize target enrichment between Plasmid library and the Selected population. Plots Extended Fig. E1j
echo -e "Rscript ${script_dir}FigE1h-i_TFi_Selected_rank.R $wd"
Rscript ${script_dir}FigE1h-i_TFi_Selected_rank.R $wd
