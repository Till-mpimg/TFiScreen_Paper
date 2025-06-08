#!/bin/bash
# Master script to analyze qPCR data and plot Extended Fig. E9a-c, h

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Processes and visualizes qPCR data following Xist activator knockdowns for Extended Fig. E9a
echo -e "Rscript ${script_dir}FigE9a_XistAct_qPCR.R $wd"
Rscript ${script_dir}FigE9a_XistAct_qPCR.R $wd

#Processes and visualizes qPCR data from the CUT&Tag/FISH experiment for Extended Fig. E9b-c
echo -e "Rscript ${script_dir}FigE9bc_qPCR_cont.R $wd"
Rscript ${script_dir}FigE9bc_qPCR_cont.R $wd

#Processes and visualizes qPCR data investigating effect of Xist activator knockdowns on differentiation marker genes Extended Fig. E9h
echo -e "Rscript ${script_dir}FigE9h_DIFFmarkers_qpcr.R $wd"
Rscript ${script_dir}FigE9h_DIFFmarkers_qpcr.R $wd



