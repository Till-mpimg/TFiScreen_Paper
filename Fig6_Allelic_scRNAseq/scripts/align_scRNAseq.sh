#!/bin/bash

###############################
## ~~~ STARsolo       ####
##
###############################
# this script runs STARsolo with allele-specific mapping, WASP filtering (improves allele-specific mapping)
# using the mm10_edit reference genome (mm10 from cellranger "refdata-gex-mm10-2020-A", plus added annotation for: Xert, Linx, CasTuner transgene, multiguide transgene N-masked protospacer)


wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'
fastq_dir=${input_dir}scRNAseq_fast'/'



## Run STARSolo 
genomeDir=${data_dir}STAR_Masked
soloCBwhitelist=${input_dir}cellbarcode_whitelist.txt
vcfFile=${input_dir}TX1072_SNPs.vcf

cd ${fastq_dir}

for f in $(ls scRNAseq_rep*_R1.fastq.gz | rev | cut -c 13- | rev | uniq)
do	
	STAR \
	--genomeDir $genomeDir \
	--readFilesIn $f\_R1.fastq.gz $f\_R2.fastq.gz \
	--readFilesCommand zcat \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist $soloCBwhitelist \
	--soloFeatures GeneFull_Ex50pAS \
	--soloUMIlen 12 \
	--soloCellFilter EmptyDrops_CR \
	--outFileNamePrefix $data_dir \
	--varVCFfile $vcfFile \
	--waspOutputMode SAMtag \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vG vA vW ha \
	--genomeTransformOutput SAM SJ Quant \
	--outSAMtype BAM SortedByCoordinate
done


#Merge BAM files together for downstream analysis
bams=$(find $data_dir -type f -name "*_Aligned.sortedByCoord.out.bam" | xargs)

samtools merge -o ${output_dir}scRNAseq_combi.bam $bams
samtools index ${output_dir}scRNAseq_combi.bam
