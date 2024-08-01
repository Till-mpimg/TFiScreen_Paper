#!/bin/bash
#Script to align RNAseq timecourse data for figure 2

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'


fastq_dir=${input_dir}rnaseq_fastq'/' #Specify fastq directory
genome_dir=${input_dir}STAR_masked'/' #Specify mm10 genome (masked for TX SNPs)

cd $data_dir


mkdir -p final_bam
bam_dir=${data_dir}/final_bam'/'

cd $fastq_dir

for f in $(ls *.fastq.gz | rev | cut -c 12- | rev | uniq)
do
	cd ${data_dir}
	
	#Trim read 1 based on quality and polyA sections
	echo -e "Trimming $f with trim_galore"
	prun python3 trim_galore --basename $f ${fastq_dir}$f\_1.fastq.gz  >> $f\.trimmingStats.txt 2>&1
	prun python3 trim_galore --polyA --basename $f\_polyA $f\_trimmed.fq.gz >> $f\_polyA.trimmingStats.txt 2>&1
	
	#Align Read 1 single-end
	echo -e "Mapping $f with STAR"
	STAR --genomeDir $genome_dir --readFilesIn $f\_polyA_trimmed_trimmed.fq.gz --readFilesCommand zcat \
				--outSAMtype BAM Unsorted --outFileNamePrefix $f\_ --outSAMattributes NH HI NM MD
	
	#Sort and index resulting BAM files
	echo -e "Sorting BAM files for $f"
	samtools sort -m 1G $f\_Aligned.out.bam -T $f\_sorted -o ${bam_dir}$f\_sorted.bam
	samtools index ${bam_dir}$f\_sorted.bam

done
