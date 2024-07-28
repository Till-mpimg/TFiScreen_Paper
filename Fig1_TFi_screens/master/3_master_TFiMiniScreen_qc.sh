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

#Processes FACS data and plots Extended Fig. E3d
echo -e "Rscript ${script_dir}FigE3d_TFiMini_FACS.R $wd"
Rscript ${script_dir}FigE3d_TFiMini_FACS.R $wd

#Plot qPCR data for Extended Fig. E3e
echo -e "Rscript ${script_dir}FigE3e_TFiMini_qPCR.R $wd"
Rscript ${script_dir}FigE3e_TFiMini_qPCR.R $wd

#Visualize QC for Extended Fig. E3f-g
echo -e "Rscript ${script_dir}FigE3f-g_TFiMini_QC.R $wd"
Rscript ${script_dir}FigE3f-g_TFiMini_QC.R $wd

#Visualize target enrichment between control populations. Plots Extended Fig. E4e-f
echo -e "Rscript ${script_dir}FigE4e-f_TFiMini_Selected_rank.R $wd"
Rscript ${script_dir}FigE4e-f_TFiMini_Selected_rank.R $wd
