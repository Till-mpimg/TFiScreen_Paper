#!/bin/bash
# Master script to design TFi library and plot Extended Fig. E1b-d

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Design TFi Library and plot Figs. E1b
echo -e "Rscript ${script_dir}FigE1b_TFi_Lib_Design.R $wd"
Rscript ${script_dir}FigE1b_TFi_Lib_Design.R $wd

#Design guides from Guidescan2 output
echo -e "Rscript ${script_dir}TFi_Guide_Design.R $wd"
Rscript ${script_dir}TFi_Guide_Design.R $wd

#Plot library design info
echo -e "Rscript ${script_dir}FigE1c-d_Plot_TFi_Design.R $wd"
Rscript ${script_dir}FigE1c-d_Plot_TFi_Design.R $wd

