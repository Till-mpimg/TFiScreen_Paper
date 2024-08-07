#!/bin/bash


# Get input file and locations  
wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'
bam=${output_dir}scRNAseq_combi.bam
genes_gtf=${input_dir}GENCODE_vM25_plus_Xert_Linx.gtf

### Non-allele-specific counting

# Remove missing cell barcodes because umi_tools can't deal with this
samtools view -h $bam | grep -v "CB:Z:-" | samtools view -b -o $data_dir'scRNAseq_combi_barcodesFiltered.bam'
samtools index $data_dir'scRNAseq_combi_barcodesFiltered.bam'

# Assign genes with featurecounts and sort
featureCounts -a $genes_gtf -o  $data_dir'scRNAseq_combi_featureCounts_output.txt' -R BAM $data_dir'scRNAseq_combi_barcodesFiltered.bam' -T 16 -s 1
samtools sort $data_dir'scRNAseq_combi_barcodesFiltered.bam.featureCounts.bam' -o  $output_dir'scRNAseq_combi_assigned_sorted.bam'
samtools index $output_dir'scRNAseq_combi_assigned_sorted.bam'

# Count with umi tools
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $output_dir'scRNAseq_combi_assigned_sorted.bam' -S $output_dir'scRNAseq_combi_counts.tsv.gz' --umi-tag=UB --cell-tag=CB --extract-umi-method=tag

# Unzip
gunzip $output_dir'scRNAseq_combi_counts.tsv.gz'



### Allele-specific counting

# Split bams haplotypes depending on ha-tag and filter for WASP
samtools view -H $output_dir'scRNAseq_combi_assigned_sorted.bam' > $data_dir'header.sam'
samtools view $output_dir'scRNAseq_combi_assigned_sorted.bam'  | grep -E '\s(vW:i:1.*ha:i:1|ha:i:1.*vW:i:1)\s' > $data_dir'scRNAseq_combi_assigned_sorted_haplo1.sam'
cat $data_dir'header.sam' $data_dir'scRNAseq_combi_assigned_sorted_haplo1.sam' | samtools view -b -o $data_dir'scRNAseq_combi_assigned_sorted_haplo1.bam'
samtools index $data_dir'scRNAseq_combi_assigned_sorted_haplo1.bam'
samtools view $data_dir'scRNAseq_combi_assigned_sorted.bam' | grep -E '\s(vW:i:1.*ha:i:2|ha:i:2.*vW:i:1)\s' > $data_dir'scRNAseq_combi_assigned_sorted_haplo2.sam'
cat $data_dir'header.sam' $data_dir'scRNAseq_combi_assigned_sorted_haplo2.sam' | samtools view -b -o $data_dir'scRNAseq_combi_assigned_sorted_haplo2.bam'
samtools index $data_dir'scRNAseq_combi_assigned_sorted_haplo2.bam'

# Extract X chromosomal reads
samtools view $data_dir'scRNAseq_combi_assigned_sorted_haplo1.bam' chrX -b > $name'_assigned_sorted_haplo1_chrX.bam'
samtools view $data_dir'scRNAseq_combi_assigned_sorted_haplo2.bam' chrX -b > $name'_assigned_sorted_haplo2_chrX.bam'
samtools index $data_dir'scRNAseq_combi_assigned_sorted_haplo1_chrX.bam'
samtools index $data_dir'scRNAseq_combi_assigned_sorted_haplo2_chrX.bam'

# Clean-up
rm $bampath'header.sam'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo1.sam'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo2.sam'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo1.bam'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo2.bam'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo1.bam.bai'
rm $data_dir'scRNAseq_combi_assigned_sorted_haplo2.bam.bai'

# Count X-chromosomal reads
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $data_dir'scRNAseq_combi_assigned_sorted_haplo1_chrX.bam' -S $output_dir'scRNAseq_combi_counts_haplo1_chrX.tsv.gz' --umi-tag=UB --cell-tag=CB --extract-umi-method=tag
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $data_dir'scRNAseq_combi_assigned_sorted_haplo2_chrX.bam' -S $output_dir'scRNAseq_combi_counts_haplo2_chrX.tsv.gz' --umi-tag=UB --cell-tag=CB --extract-umi-method=tag

# Unzip
gunzip $output_dir'scRNAseq_combi_counts_haplo1_chrX.tsv.gz'
gunzip $output_dir'scRNAseq_combi_counts_haplo2_chrX.tsv.gz'








