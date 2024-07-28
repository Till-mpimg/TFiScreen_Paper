#!/bin/bash
#Runs mageck commands for the full TFi Lib (also corrects the duplicated guides)

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'

cd ${input_dir}

files=$(ls -p *TFi*_1.fastq.gz | grep -v / | tr '\n' ' ')
names=$(echo -e "$files" | sed "s/XX_D2_TFi_//g" | sed "s/_1.fastq.gz//g" | sed "s/ /,/g" |  sed "s/.$//g")

echo -e "Counting Reads via Mageck"
mageck count -l TFi_library_file.txt --fastq $files \
  --sample-label $names \
  --norm-method control --control-sgrna TFi_NT_controls.txt --pdf-report -n ${data_dir}TFi

mv ${data_dir}TFi_countsummary.pdf ${output_dir}TFi_countsummary.pdf

echo -e "Correcting duplicate guides"
Rscript ${wd}scripts/TFi_Correct_Guides.R $data_dir $input_dir $output_dir

echo -e "Running Mageck MLE"
echo -e "Comparing High and Negative fractions"
mageck mle -k ${output_dir}TFi.count_corr_normalized.txt -d ${input_dir}TFi_design_matrix_HvN.txt -n ${data_dir}TFi_mle_HvN \
  --permutation-round 10

mv ${data_dir}TFi_mle_HvN.gene_summary.txt ${output_dir}TFi_mle_HvN.gene_summary.txt

echo -e "Comparing Low and Negative fractions"
mageck mle -k ${output_dir}TFi.count_corr_normalized.txt -d ${input_dir}TFi_design_matrix_LvN.txt -n ${data_dir}TFi_mle_LvN \
  --permutation-round 10

mv ${data_dir}TFi_mle_LvN.gene_summary.txt ${output_dir}TFi_mle_LvN.gene_summary.txt

echo -e "Comparing Selection fractions and Plasmid Library"
mageck mle -k ${output_dir}TFi.count_corr_normalized.txt -d ${input_dir}TFi_design_matrix_SvP.txt -n ${data_dir}TFi_mle_SvP \
  --permutation-round 10

mv ${data_dir}TFi_mle_SvP.gene_summary.txt ${output_dir}TFi_mle_SvP.gene_summary.txt
