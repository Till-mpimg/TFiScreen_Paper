#!/bin/bash
# Master script to design TFiMini library, generate count tables, enrich screen targets and perform QC analysis. Plot Extended Fig. E3d-g, E4e-f.

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Design TFiMini Library
echo -e "Rscript ${script_dir}TFiMini_Lib_Design.R $wd"
Rscript ${script_dir}TFiMini_Lib_Design.R $wd

#Performs counting and target enrichment using MAGeCK
echo -e "${script_dir}TFiMini_mageck.sh $wd"
${script_dir}TFiMini_mageck.sh $wd

#Processes FACS data and plots Extended Fig. E3a
echo -e "Rscript ${script_dir}FigE3a_TFiMini_FACS.R $wd"
Rscript ${script_dir}FigE3a_TFiMini_FACS.R $wd

#Calculates fold change between sorting gates for TFi and TFiMini screens and plots Extended Fig. E3b
echo -e "Rscript ${script_dir}FigE3b_Screen_sorting_width.R $wd"
Rscript ${script_dir}FigE3b_Screen_sorting_width.R $wd

#Plot qPCR data for Extended Fig. E3c
echo -e "Rscript ${script_dir}FigE3c_TFiMini_qPCR.R $wd"
Rscript ${script_dir}FigE3c_TFiMini_qPCR.R $wd

#Visualize QC for Extended Fig. E3d-e
echo -e "Rscript ${script_dir}FigE3d-e_TFiMini_QC.R $wd"
Rscript ${script_dir}FigE3d-e_TFiMini_QC.R $wd

#Visualize target enrichment between control populations. Plots Extended Fig. E3f-g
echo -e "Rscript ${script_dir}FigE3f-g_TFiMini_Selected_rank.R $wd"
Rscript ${script_dir}FigE3f-g_TFiMini_Selected_rank.R $wd
