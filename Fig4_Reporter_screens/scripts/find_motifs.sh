#!/bin/bash
#Finds motifs in REs from TFiMini library

wd=$1
data_dir=${wd}data'/'
input_dir=${wd}input_files'/'
output_dir=${wd}output_files'/'

meme=${wd}input_files/Fig_E7/consensus_pwms_jvierstra.meme
TFiMini_genes=${wd}input_files/Fig_E7/TFiMini_genes.txt
re_total=${wd}input_files/Fig_E7/re_total.bed
mm10=${wd}input_files/Fig_E7/mm10.fa
jvierstra=${wd}input_files/Fig_E7/jvierstra_motif_metadata.tsv

cd $data_dir

#Writes fasta file to search from the mm10 genome and the RE coordinates
echo -e "getfasta for re_total"
bedtools getfasta -fi $mm10 -bed $re_total > ${data}re_total.fa

echo -e "fimo for motif set"
fimo -oc ${output_dir} $meme ${data}re_total.fa

#Gives out all motifs in the REs
echo -e "Creating list of interesting motifs"
> ${output_dir}TFiMini_motifs.txt
TFiMini_motifs=${output_dir}TFiMini_motifs.txt

#Use this code for jVierstra motifs
#Translates human motifs (ZNF) to mouse factors (ZFP)
while read gene
do
  if [[ $gene =~ [Zz][Ff][Pp] ]]; then
    number=$(grep -oPi 'ZFP\K\d+' <<< "$gene")
    pattern="(ZFP$number|ZNF$number)"
    grep -Pi "\t$pattern\t" $jvierstra >> $TFiMini_motifs
  else
    grep -Pi "\t$gene\t" $jvierstra >> $TFiMini_motifs
  fi
done < $TFiMini_genes
