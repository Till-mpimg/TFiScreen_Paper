#!/bin/bash
# Puts all count tables together and creates metasheet

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'


output=${data_dir}mohammed_2017_counts.txt
names=${input_dir}Mohammed_cell_names.txt

mkdir -p ${data_dir}Mohammed_single_cells'\'
count_dir=${data_dir}Mohammed_single_cells'\'

while IFS=$'\t' read -r sra name; do
    cp ${data_dir}$sra\_ReadsPerGene.out.tab ${count_dir}$name\_counts.txt
done < "$names"


echo -e "GENCODE" > $output
sed 1,4d ${count_dir}E3.5_P8_Cell2_embryo1_single_counts.txt | cut -f1 -d$'\t' - >> ${data_dir}mohammed_2017_counts.txt

for f in $(ls *_counts.txt | rev | cut -c 12- | rev | uniq)
do
	echo -e "$f" > ${data_dir}$f\_temp_col.txt
	sed 1,4d $f\_counts.txt | cut -f2 -d$'\t' - >> ${data_dir}$f\_temp_col.txt
	paste $output ${data_dir}$f\_temp_col.txt > ${data_dir}$f\_temp.txt
	mv ${data_dir}$f\_temp.txt $output
done

