#!/bin/bash
# Master script to process Reporter screen data and plot quality control metrics for Fig. E8a-e, E9a-b

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Analyzes FACS data from the reporter screens and plot Extended Fig. E8a-b
echo -e "Rscript ${script_dir}FigE8ab_reporter_screen_facs.R $wd"
Rscript ${script_dir}FigE8ab_reporter_screen_facs.R $wd

#Counts guides in the reporter screens using mageck
echo -e "${script_dir}Count_reporter_screens.sh $wd"
${script_dir}Count_reporter_screens.sh $wd

#Perform quality control on the counts
echo -e "Rscript ${script_dir}FigE8cde_REPi_QC.R $wd"
Rscript ${script_dir}FigE8cde_REPi_QC.R $wd

#Quantify noRE screen, calculate interaction scores and plot Fig. 4f-g and Extended Fig. E9b
echo -e "Rscript ${script_dir}Fig4fgE9a_Calc_interaction_score.R $wd"
Rscript ${script_dir}Fig4fgE9a_Calc_interaction_score.R $wd

#Quantify enrichment of positive controls and plot Fig. E9a
echo -e "Rscript ${script_dir}FigE9a_REPi_pos_cont.R $wd"
Rscript ${script_dir}FigE9a_REPi_pos_cont.R $wd



