#!/bin/bash
# Master script to process scRNAseq data and produce count tables

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Analyze qPCR data and plot Extended Fig. E10b
echo -e "Rscript ${script_dir}FigE10b_qPCR_CRISPRi_cont.R $wd"
Rscript ${script_dir}FigE10b_qPCR_CRISPRi_cont.R $wd

#Segments and quantifies FISH data for the sgFtx-Xert-REs control experiment
echo -e "prun python3 ${script_dir}CRISPRi_cont_FISH_segmentation.py $wd"
prun python3 ${script_dir}CRISPRi_cont_FISH_segmentation.py $wd

#Plots FISH data for Extended Fig. E10c-d
echo -e "Rscript ${script_dir}FigE10cd_FISH_CRISPRi_cont.R $wd"
Rscript ${script_dir}FigE10cd_FISH_CRISPRi_cont.R $wd



