#!/bin/bash
#Merges CUT&Tag tracks for visualization

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
out_dir=${wd}output_files'/'

cd $data_dir

for f in $(ls *_R1_dedup.bam | rev | cut -c 14- | rev | uniq)
do
	echo -e "Processing $f"
	
	samtools merge -f ${out_dir}${f}.bam ${f}_R1_dedup.bam ${f}_R2_dedup.bam ${f}_R3_dedup.bam 
	samtools index ${out_dir}${f}.bam
	
	prun python3 bamCoverage -b ${out_dir}${f}.bam  -o ${out_dir}${f}.bw -bs 10 -e --normalizeUsing CPM -ignore chrX chrY
done
