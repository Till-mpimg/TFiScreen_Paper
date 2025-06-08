#!/bin/bash
# Master script to process scRNAseq data and produce count tables

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

mkdir -p ${wd}data'/'
mkdir -p ${wd}output_files'/'

#Build N-masked TX1072 genome for STAR solo (GENCODE supplemented with Xert/Linx coordinates)
echo -e "${script_dir}index_genome_starsolo.sh $wd"
${script_dir}index_genome_starsolo.sh $wd

#Aligns fastq files (subset for matching cell barcodes/sample) to mm10 with STAR solo
echo -e "${script_dir}align_scRNAseq.sh $wd"
${script_dir}align_scRNAseq.sh $wd

#Creates total count tables
echo -e "${script_dir}create_count_tables.sh $wd"
${script_dir}create_count_tables.sh $wd

#Creates total count tables
echo -e "prun python3 ${script_dir}Fig6bE11abc_prepare_scData.py $wd"
prun python3 ${script_dir}Fig6bE11abc_prepare_scData.py $wd"



