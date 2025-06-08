#!/bin/bash

# Runs STARsolo genome indexing
# with the "diploid genome" function that incorporates genomic variants in the reference - useful for allele-specific XCI analysis
# Constructs annotation with Xert, Linx, N-masked SNPs and CasTuner construct

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'

# Before aligning your reads, you need to generate an index of your reference genome.
# 100 bp read length, RAM 128 gb = 137438953472 bytes
STAR \
--runMode genomeGenerate \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM 25000000000 \
--genomeTransformVCF ${input_dir}TX1072_SNPs.vcf \
--genomeTransformType Diploid \
--genomeSuffixLengthMax 101 \
--genomeDir "${data_dir}STAR_Masked" \
--genomeFastaFiles "${input_dir}mm10_edit.fa" \
--sjdbGTFfile "${input_files}R/GENCODE_edit.gtf"
