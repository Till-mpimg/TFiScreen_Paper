#!/bin/bash
# Integrates in vitro RNA-timecourse data with in vivo scRNAseq (from Mohammed 2017)

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'


#Retrieves fastq files from Mohammed 2017
echo -e "${script_dir}Mohammed_fastq_dump.sh $wd"
${script_dir}Mohammed_fastq_dump.sh $wd

#Aligns Muhammed 2017 single cell data
echo -e "${script_dir}Align_Mohammed2017.sh $wd"
${script_dir}Align_Mohammed2017.sh $wd

#Combines count tables into count matrix
echo -e "${script_dir}Mohammed_generate_count_matrix.sh $wd"
${script_dir}Mohammed_generate_count_matrix.sh $wd

#Performs DESeq2 comparing 0h with other timepoints
echo -e "Rscript ${script_dir}perform_DESeq2_timescale.R $wd"
Rscript ${script_dir}perform_DESeq2_timescale.R $wd

#Does the integration with the scATAcat pipeline
echo -e "prun python3 ${script_dir}Integrate_invivo.R $wd"
prun python3 ${script_dir}Integrate_invivo.py $wd

#Plots the PCA results for Fig. 3b 
echo -e "Rscript ${script_dir}Fig3b_plot_invivo_pca.R $wd"
Rscript ${script_dir}Fig3b_plot_invivo_pca.R $wd

#Plots the distance matrix for Extended Fig. E5a
echo -e "Rscript ${script_dir}FigE5a_plot_distance_matrix.R $wd"
Rscript ${script_dir}FigE5a_plot_distance_matrix.R $wd

