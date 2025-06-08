#!/bin/bash
# Master script to process CUT&Tag data, produce bigwig tracks and plot Fig. 5c and Extended Fig. E9e-g

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Processes and aligns CUT&Tag data
echo -e "${script_dir}CUTnTag_align.sh $wd"
${script_dir}CUTnTag_align.sh $wd

#Merges replicate CUT&Tag BIGWIG tracks
echo -e "${script_dir}CUTnTag_merge.sh $wd"
${script_dir}CUTnTag_merge.sh $wd

#Plots enrichment of CUT&Tag data at reporter screen REs
echo -e "Rscript ${script_dir}Fig5c_RE_CnT_XistAct.R $wd"
Rscript ${script_dir}Fig5c_RE_CnT_XistAct.R $wd

#Plots comparison between published Chip data and reporter screens for Extended Fig. E9e
echo -e "Rscript ${script_dir}FigE9e_chip_reporter_comp.R $wd"
Rscript ${script_dir}FigE9e_chip_reporter_comp.R $wd

#Call differential peaks between sgNT and Xist activator knockdowns
echo -e "Rscript ${script_dir}DiffBind_CnT_XistAct.R $wd"
Rscript ${script_dir}DiffBind_CnT_XistAct.R $wd

#Performs genome-wide CUT&Tag analysis and plots Extended Fig. E9f
echo -e "Rscript ${script_dir}FigE9f_CnT_scatter.R $wd"
Rscript ${script_dir}FigE9f_CnT_scatter.R $wd

#Plots euler diagrams for Extended Fig. E9g
echo -e "Rscript ${script_dir}FigE9g_CnT_euler.R $wd"
Rscript ${script_dir}FigE9g_CnT_euler.R $wd
