#!/bin/bash
#Script to retrieve FASTQ files for Mohammed 2017

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'

mkdir -p ${data_dir}Mohammed_fastq'/'
fastq_dir=${data_dir}Mohammed_fastq'/'

acc_list=${input_dir}Mohammed_SRR_Acc_List.txt

cd ${fastq_dir}
while read line
do
  prefetch $line && fasterq-dump $line
done < $acc_list
