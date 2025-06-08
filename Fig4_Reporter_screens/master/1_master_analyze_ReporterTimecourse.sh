#!/bin/bash
# Master script to process FACS data from reporter time course and plot Fig. 4c and Extended Fig. E6b-e

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Processes and visualizes reporter time course FACS data
echo -e "Rscript ${script_dir}Fig4cE6bc_Firewach.R $wd"
Rscript ${script_dir}Fig4cE6bc_Firewach.R $wd

#Compares reporter data with endogenous ATAC and H3K27ac
echo -e "Rscript ${script_dir}FigE6de_CnT_Comp.R $wd"
Rscript ${script_dir}FigE6de_CnT_Comp.R $wd



