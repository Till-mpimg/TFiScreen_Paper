#!/bin/bash
# Master script to analyze qPCR data and plot Extended Fig. E10a-c

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Processes and visualizes qPCR data following Xist activator knockdowns for Extended Fig. E10a
echo -e "Rscript ${script_dir}FigE10a_XistAct_qPCR.R $wd"
Rscript ${script_dir}FigE10a_XistAct_qPCR.R $wd

#Processes and visualizes qPCR data from the CUT&Tag/FISH experiment for Extended Fig. E10b-c
echo -e "Rscript ${script_dir}FigE10bc_qPCR_cont.R $wd"
Rscript ${script_dir}FigE10bc_qPCR_cont.R $wd



