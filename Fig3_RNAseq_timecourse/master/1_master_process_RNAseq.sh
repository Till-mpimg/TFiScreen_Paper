#!/bin/bash
# Master script to process RNA-seq time course data, quantify gene expression and group TF genes by expression dynamics.

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'


#Aligns RNA-seq data
echo -e "${script_dir}RNA_Timecourse_align.sh $wd"
${script_dir}RNA_Timecourse_align.sh $wd

#Uses Rsubread to quantify read counts and calculates CPM
echo -e "Rscript ${script_dir}Timecourse_CPM_calc.R $wd"
Rscript ${script_dir}Timecourse_CPM_calc.R $wd

#Performs DESeq2 analysis per timepoint
echo -e "Rscript ${script_dir}perform_DESeq.R $wd"
Rscript ${script_dir}perform_DESeq.R $wd

#Groups TFi genes by expression dynamics
echo -e "Rscript ${script_dir}cluster_TF_timecourse.R $wd"
Rscript ${script_dir}cluster_TF_timecourse.R $wd



