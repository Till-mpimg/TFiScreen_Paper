#!/bin/bash
#Merges CUT&Tag tracks for visualization
#Merges peak files for downstream analysis

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
out_dir=${wd}output_files'/'

cd $data_dir

#Merging .bam and .narrowPeak files by samples
for f in $(ls *_R1_dedup.bam | rev | cut -c 14- | rev | uniq)
do
	echo -e "Processing $f"
	
	samtools merge -f ${out_dir}${f}.bam ${f}_R1_dedup.bam ${f}_R2_dedup.bam ${f}_R3_dedup.bam 
	samtools index ${out_dir}${f}.bam
	
	prun python3 bamCoverage -b ${out_dir}${f}.bam  -o ${out_dir}${f}.bw -bs 10 -e --normalizeUsing CPM -ignore chrX chrY
	
	bedtools intersect -a ${out_dir}${f}-R1_dedup_peaks.narrowPeak -b ${out_dir}${f}-R2_dedup_peaks.narrowPeak \
      > ${data_dir}${f}_merged_R12.narrowPeak
      
    bedtools intersect -a ${data_dir}${f}_merged_R12.narrowPeak -b ${out_dir}${f}-R3_dedup_peaks.narrowPeak \
      > ${data_dir}${f}_merged.narrowPeak
      
	awk 'BEGIN {OFS="\t"}; {print $1,$2,$3}' ${data_dir}${f}_merged.narrowPeak \
      > ${out_dir}${f}.bed
      
    rm ${data_dir}${f}_merged.narrowPeak
    rm ${data_dir}${f}_merged_R12.narrowPeak
done

#Creating composite peaks files
cd $out_dir

cat XX_D2_H3K27ac_sgNT.bed XX_D2_H3K27ac_sgOtx2.bed XX_D2_H3K27ac_sgZic3.bed XX_D2_H3K27ac_sgNfrkb.bed XX_D2_H3K27ac_sgFoxd3.bed | \
	bedtools sort -i - | bedtools merge -i - > ${out_dir}H3K27ac_composite.bed

cat XX_D2_H3K4me3_sgNT.bed XX_D2_H3K4me3_sgOtx2.bed XX_D2_H3K4me3_sgZic3.bed H3K4me3-sgNfrkb.bed XX_D2_H3K4me3_sgFoxd3.bed | \
	bedtools sort -i - | bedtools merge -i - > ${out_dir}H3K4me3_composite.bed

cat XX_D2_H3K4me1_sgNT.bed XX_D2_H3K4me1_sgOtx2.bed XX_D2_H3K4me1_sgZic3.bed XX_D2_H3K4me1_sgNfrkb.bed XX_D2_H3K4me1_sgFoxd3.bed | \
	bedtools sort -i - | bedtools merge -i - > ${out_dir}H3K4me1_composite.bed
