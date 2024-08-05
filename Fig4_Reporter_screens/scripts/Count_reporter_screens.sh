#!/bin/bash
#Runs mageck commands for the full REPi Lib 

wd=$1

data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'
fastq_dir=${input_dir}REPi_fastq/

cd $fastq_dir

files=$(ls -p *_1.fastq.gz | grep -v / | tr '\n' ' ')
names=$(echo -e "$files" | sed "s/XX_D2_REPi_//g" | sed "s/_1.fastq.gz//g"| sed "s/GFP-//g"  | sed "s/ /,/g" |  sed "s/.$//g")

#Counts sgRNA occurences in read 1 and normalizes by the NT controls
mageck count -l ${input_dir}REPi_library_file.txt --fastq $files --sample-label $names --norm-method control --control-sgrna ${input_dir}REPi_NT_file.txt  --pdf-report -n ${data_dir}REPi
mv ${data_dir}REPi_countsummary.pdf ${output_dir}REPi_countsummary.pdf
mv ${data_dir}REPi.countsummary.txt ${output_dir}REPi.countsummary.txt
mv ${data_dir}REPi.count.txt ${output_dir}REPi.count.txt
mv ${data_dir}REPi.count_normalized.txt ${output_dir}REPi.count_normalized.txt
