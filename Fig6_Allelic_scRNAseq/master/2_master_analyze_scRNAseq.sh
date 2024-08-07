#!/bin/bash
# Master script to analyze scRNAseq for chrX allelic expression

dir=$(pwd)
wd="$(dirname "$dir")"'/'
script_dir=${wd}scripts'/'

#Filters for monoallelic Xist-positive and plots Fig. 6d and Extended Fig. E11d
echo -e "Rscript ${script_dir}Fig6cdE11d_generate_filtered_table.R $wd"
Rscript ${script_dir}Fig6dE11d_generate_filtered_table.R $wd

#Analyzes XCI in groups of cells with variable Xist levels. Plots Fig. 6f-h and Extended Fig. E11f-i
echo -e "Rscript ${script_dir}Fig6fghE11fghi_wt_xist_bins.R $wd"
Rscript ${script_dir}Fig6fghE11fghi_wt_xist_bins.R $wd

#Compares XCI between deletion and wildtype cells. Plots Fig. 6e and Extended Fig. E11e
echo -e "${script_dir}create_count_tables.sh $wd"
${script_dir}create_count_tables.sh $wd



