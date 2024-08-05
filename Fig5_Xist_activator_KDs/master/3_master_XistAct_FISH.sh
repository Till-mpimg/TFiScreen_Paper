#!/bin/bash
# Master script to process FISH data and plot Fig. 5-g and Extended Fig. E10d

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Segments FISH data and produces output tables (intensity/frequency of Xist clouds). Also produces images for Fig. 5d
echo -e "prun python3 ${script_dir}Fig5d_FISH_segmentation.py $wd"
prun python3 ${script_dir}Fig5d_FISH_segmentation.py $wd

#Plots FISH output as Fig. 5e-g and Extended Fig. E10d 
echo -e "Rscript ${script_dir}Fig5efgE10d_plot_FISH.R $wd"
Rscript ${script_dir}Fig5efgE10d_plot_FISH.R $wd



