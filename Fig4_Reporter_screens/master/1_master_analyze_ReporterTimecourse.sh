#!/bin/bash
# Master script to process FACS data from reporter time course and plot Fig. 4c and Extended Fig. E7b-e

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Processes and visualizes reporter time course FACS data
echo -e "Rscript ${script_dir}Fig4cE7bc_Firewach.R $wd"
Rscript ${script_dir}Timecourse_CPM_calc.R $wd

#Compares reporter data with endogenous ATAC and H3K27ac
echo -e "Rscript ${script_dir}FigE7de_CnT_Comp.R $wd"
Rscript ${script_dir}FigE7de_CnT_Comp.R $wd



