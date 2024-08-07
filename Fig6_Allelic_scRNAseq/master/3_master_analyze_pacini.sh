#!/bin/bash
# Master script to analyze pacini scRNAseq for chrX allelic expression

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Filters Xist positive Day 3 cells from Pacini et al. 2021
echo -e "Prun python3 ${script_dir}analyze_pacini.py $wd"
Rscript ${script_dir}Fig6dE11d_generate_filtered_table.R $wd

#Repeats Xist bin analysis for Pacini data and plots Extended Fig. E11j
echo -e "Rscript ${script_dir}FigE11j_Xist_bins_pacini.R $wd"
Rscript ${script_dir}FigE11j_Xist_bins_pacini.R $wd




