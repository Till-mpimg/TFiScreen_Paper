#!/bin/bash
# Master script to analyze the TFi screens and plot Fig. 1b-c and Extended Fig. E2a-c, E4a-d

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Plots rank plot for High/Negative comp in the TFi screen as Fig. 1b
echo -e "Rscript ${script_dir}Fig1b_TFi_High_rank.R $wd"
Rscript ${script_dir}Fig1b_TFi_High_rank.R $wd

#Plots beta score heatmaps for TFi High/Negative, TFi Low/Negative and TFiMini High/Negative as Fig. 1c
echo -e "Rscript ${script_dir}Fig1b_TFi_Low_High_heat.R $wd"
Rscript ${script_dir}Fig1b_TFi_Low_High_heat.R $wd

#Plots results vs Unsorted for the TFi lib as Extended Fig. E2a
echo -e "Rscript ${script_dir}FigE2a_TFi_vU_heatmap.R $wd"
Rscript ${script_dir}FigE2a_TFi_vU_heatmap.R $wd

#Plots per guide results for the TFi lib as Extended Fig. E2b
echo -e "Rscript ${script_dir}FigE2b_TFi_guide_plots.R $wd"
Rscript ${script_dir}FigE2b_TFi_guide_plots.R $wd

#Plots rank plot for Low/Negative comp in the TFi screen as Extended Fig. E2c
echo -e "Rscript ${script_dir}FigE2c_TFi_Low_rank.R $wd"
Rscript ${script_dir}FigE2c_TFi_Low_rank.R $wd

#Plots rank plot for High/Negative comp in the TFiMini screen as Extended Fig. E4a
echo -e "Rscript ${script_dir}FigE4a_TFiMini_High_rank.R $wd"
Rscript ${script_dir}FigE4a_TFiMini_High_rank.R $wd

#Plots scatter plot and venn diagram to compare the TFi and TFiMini screens as Extended Fig. E4b-d
echo -e "Rscript ${script_dir}FigE4b-d_screens_comp.R $wd"
Rscript ${script_dir}FigE4b-d_screens_comp.R $wd

