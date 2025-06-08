#!/bin/bash
# Master script to process Reporter screen data and plot quality control metrics for Fig. E8a-e, E9a-b

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'


#Analyzes FACS data from the reporter screens and plot Extended Fig. E7a-b
echo -e "Rscript ${script_dir}FigE7ab_reporter_screen_facs.R $wd"
Rscript ${script_dir}FigE7ab_reporter_screen_facs.R $wd

#Counts guides in the reporter screens using mageck
echo -e "${script_dir}Count_reporter_screens.sh $wd"
${script_dir}Count_reporter_screens.sh $wd

#Perform quality control on the counts
echo -e "Rscript ${script_dir}FigE7cde_REPi_QC.R $wd"
Rscript ${script_dir}FigE7cde_REPi_QC.R $wd

#Quantify noRE screen, calculate interaction scores and plot Fig. 4e/g and Extended Fig. E8b
echo -e "Rscript ${script_dir}Fig4egE8a_Calc_interaction_score.R $wd"
Rscript ${script_dir}Fig4egE8a_Calc_interaction_score.R $wd

#Quantify enrichment of positive controls and plot Fig. E8a
echo -e "Rscript ${script_dir}FigE8a_REPi_pos_cont.R $wd"
Rscript ${script_dir}FigE8a_REPi_pos_cont.R $wd



