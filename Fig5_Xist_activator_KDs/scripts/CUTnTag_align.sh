#Aligns CUT&Tag of Xist activator Knockdown
#!/bin/bash

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
fastq_dir=${input_dir}CnT_fastq'/'
ebwt=${input_dir}bowtie2_mm10/mm10
blacklistmm10=${input_dir}mm10-blacklist.v2.bed
out_dir=${wd}output_files'/'

cd ${fastq_dir}

for f in $(ls *.fastq.gz | rev | cut -c 12- | rev | uniq)
do
	cd ${data_dir}
	
	echo -e "Trimming $base with trim_galore"
	prun python3 trim_galore --paired --nextera ${fastq_dir}$f\_1.fastq.gz ${fastq_dir}$f\_2.fastq.gz >> $f\.trimmingStats.txt 2>&1

	echo -e "Mapping $f with bowtie2"
	bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 2000 -p 20 -x $ebwt -1 $f\_1_val_1.fq.gz -2 $f\_2_val_2.fq.gz -S $f\.sam >> $f\.mappingStats.txt 2>&1
	rm $f\_1_val_1.fq.gz
	rm $f\_2_val_2.fq.gz

	echo -e "Creating BAM files for $f and filtering for paired and mapped reads"
	samtools view -b -h -f 2 -q 20 $f\.sam > $f\.bam
	rm $f\.sam

	echo -e "Sorting BAM files for $f"
	samtools sort -m 1G $f\.bam -T $f\_sorted -o $f\_sorted.bam

	echo -e "Removing Blacklisted regions from $f"
	bedtools intersect -v -a $f\_sorted.bam -b $blacklistmm10 > $f\_sorted_blacklisted.bam

	echo -e "Removing duplicates from $base using PICARD"
	java -jar $PICARD MarkDuplicates INPUT=$f\_sorted_blacklisted.bam OUTPUT=$f\_dedup.bam METRICS_FILE=$f\_dedupMetric.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
	samtools index $f\_dedup.bam

	echo -e "Calls peaks for $f using macs2"
	macs2 callpeak -t $f\_dedup.bam -f BAMPE -g mm -q 0.05 -n $f\ --outdir ${out_dir}

	echo -e "Generates bigwig for $base"
	prun python3 bamCoverage -b $f\_dedup.bam  -o ${out_dir}$f\.bw -bs 10 -e --normalizeUsing CPM -ignore chrX chrY
done

cd ${data_dir}
#Generating qc_file for the data
echo -e "Sample\tTotal\tFiltered\tMapped\tDuplicate\tFRiP\tPeaks" > ${out_dir}cnt_xist_act_qc_metrics.txt

for f in $(ls *_sorted_blacklisted.bam | rev | cut -c 24- | rev | uniq)
do
  echo -e "Returning qc info for $f"
  fragments=$(grep -Po -m 1 'Total reads processed.*' $f\.trimmingStats.txt | grep -Po '[0-9,]*' | tr -d ,)
  map=$(grep -Po '[0-9, /.]*% overall alignment rate' $f\.mappingStats.txt| grep -Po '[0-9,/.]*')
  mapFrac=$(printf "%.4f" $(echo $map / 100 | bc -l))
  dedup=$(printf "%.4f" $(grep -Po 'Unknown Library.*' $f\_dedupMetric.txt | awk '{print $10}'))
  filtered=$(expr $(samtools view -c $f\_dedup.bam) / 2)
  fripReads=$(expr $(bedtools intersect -a $f\_dedup.bam -b ${out_dir}$f\_peaks.narrowPeak -u -ubam | samtools view -c) / 2)
  FRiP=$(printf "%.4f" $(echo $fripReads / $filtered | bc -l))
  peaks=$(expr $(wc -l <"${out_dir}${f}_peaks.narrowPeak") - 1)

  echo -e "$f\t$fragments\t$filtered\t$mapFrac\t$dedup\t$FRiP\t$peaks" >> ${out_dir}cnt_xist_act_qc_metrics.txt
done
