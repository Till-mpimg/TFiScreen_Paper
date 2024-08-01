#!/bin/bash
# Performs alignment and filtering of scRNA-seq data

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'

ebwt=${input_dir}/STAR_masked'/'
fastq_dir=${data_dir}Mohammed_fastq'/'
gencode=${input_dir}GENCODE_vM25_plus_Xert.gtf

cd $fastq_dir

for f in $(ls *.fastq | rev | cut -c 9- | rev | uniq)
do
	cd ${data_dir}

  echo -e "Mapping $f with STAR"
  STAR --genomeDir $ebwt --readFilesIn ${fastq_dir}$f\_1.fastq ${fastq_dir}$f\_2.fastq --sjdbGTFfile $gencode \
				--outSAMtype BAM Unsorted --outFileNamePrefix $f\_ --outSAMattributes NH HI NM MD --quantMode GeneCounts
done
