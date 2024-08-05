#!/bin/bash
# Master script to process CUT&Tag data, produce bigwig tracks and plot Fig. 5c

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Processes and aligns CUT&Tag data
echo -e "${script_dir}CUTnTag_align.sh $wd"
${script_dir}CUTnTag_align.sh $wd

#Merges replicate CUT&Tag BIGWIG tracks
echo -e "${script_dir}CUTnTag_merge.sh $wd"
${script_dir}CUTnTag_merge.sh $wd

#Plots enrichment of CUT&Tag data at reporter screen REs
echo -e "Rscript ${script_dir}Fig5c_RE_CnT_XistAct.R $wd"
Rscript ${script_dir}Fig5c_RE_CnT_XistAct.R $wd



