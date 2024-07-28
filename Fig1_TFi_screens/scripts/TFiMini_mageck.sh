#!/bin/bash
#Runs mageck commands for the TFiMini Lib


wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'

cd ${input_dir}

files=$(ls -p *TFiMini*fastq.gz | grep -v / | tr '\n' ' ')
names=$(echo -e "$files" | sed "s/TFiMini_//g" | sed "s/.fastq.gz//g" | sed "s/ /,/g" |  sed "s/.$//g")



echo -e "Counting Reads via Mageck"
mageck count -l TFiMini_library_file.txt --fastq $files \
  --sample-label $names \
  --norm-method control --control-sgrna TFiMini_NT_controls.txt --pdf-report -n ${data_dir}TFiMini

mv ${data_dir}TFiMini_countsummary.pdf ${output_dir}TFiMini_countsummary.pdf
mv ${data_dir}TFiMini.count_normalized.txt ${output_dir}TFiMini.count_normalized.txt
mv ${data_dir}TFiMini.count.txt ${output_dir}TFiMini.count.txt

#echo -e "Running Mageck MLE"
echo -e "Comparing High and Negative fractions"
mageck mle -k ${output_dir}TFiMini.count_normalized.txt -d ${input_dir}TFiMini_design_matrix_HvN.txt -n ${data_dir}TFiMini_mle_HvN \
  --permutation-round 10 --max-sgrnapergene-permutation 50

mv ${data_dir}TFiMini_mle_HvN.gene_summary.txt ${output_dir}TFiMini_mle_HvN.gene_summary.txt

echo -e "Comparing Selected fractions and 2iSL Fractions"
mageck mle -k ${output_dir}TFiMini.count_normalized.txt -d ${input_dir}TFiMini_design_matrix_Sv2iSL.txt -n ${data_dir}TFiMini_mle_Sv2iSL \
  --permutation-round 10 --max-sgrnapergene-permutation 50

mv ${data_dir}TFiMini_mle_Sv2iSL.gene_summary.txt ${output_dir}TFiMini_mle_Sv2iSL.gene_summary.txt

echo -e "Comparing Selected and Plasmid Library fractions"
mageck mle -k ${output_dir}TFiMini.count_normalized.txt -d ${input_dir}TFiMini_design_matrix_SvP.txt -n ${data_dir}TFiMini_mle_SvP \
  --permutation-round 10 --max-sgrnapergene-permutation 50

mv ${data_dir}TFiMini_mle_SvP.gene_summary.txt ${output_dir}TFiMini_mle_SvP.gene_summary.txt
