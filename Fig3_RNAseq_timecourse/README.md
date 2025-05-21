# Fig. 3: Different groups of Xist activators respond to X-chromosome number and differentiation cues.

## Description
This folder contains all code necessary to reproduce the analysis for Figure 3 and Extended Figures E6. The master script has to be started from the "./master" directory. The file structure of this github should be kept intact. Any necessary files to run the master script are provided in "./input_files/" or are detailed below.


## Software dependencies and operating systems
In order to perform the analyses the following software has to be installed and available on the command line ($PATH):
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- SRA Toolkit (v3.0.7) collection of C scripts from "https://github.com/ncbi/sra-tools/"
- STAR (v2.5.4a) aligner from "https://github.com/alexdobin/STAR"
- TrimGalore (v0.6.4) PERL script from "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"

Python (v3.10.10) script can be run with the following libraries:
- anndata (v0.10.3)
- matplotlib (v3.8.2)
- numpy (v1.26.2)
- pandas (v2.1.4)
- pickle (v0.7.5)
- scanpy (v1.9.4)
- scATAcat (v0.0.1)
- scipy (v1.11.4)
- seaborn (v0.13.0) 
- sklearn (v1.3.2)

R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- DESeq2 (v1.36.0)
- egg (v0.4.5)
- EnvStats (v2.8.1)
- gg4hx (v0.2.8)
- gridExtra (v2.3)
- pheatmap (v1.0.12)
- psych (v2.4.3)
- Rgb (v1.7.5)
- Rsubread (v2.12.3)
- scales (v1.3.0)
- tidyverse (v1.3.2)
- viridis (v0.6.2)


## Reproduce analysis
All master scripts can be run independently from each other. Detailed instructions are listed below.
- "1_master_process_RNAseq.sh": Aligns RNA-seq time course data to the mm10 genome. Quantifies CPM and calls DE genes between XX and XO cells. Groups TF genes according to expression patterns. FASTQ.GZ files from the RNA-seq should be retrieved from GSE273071(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273071), named according to the pattern "RNA_[XX|XO]_[0|10|16|24|30|36|48|56|72|96]h_R[1|2|3]_[1|2].fastq.gz" and stored in "./input_files/rnaseq_fastq/". STAR genome built with TX1072 N-masked SNPs should be downloaded from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609), unzipped and stored in "./input_files/". "GENCODE_vM25_plus_Xert.gtf" should be downloaded from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".
- "2_master_integrate_invivo.sh": Performs integration of the in vitro RNA-seq timecourse with in vivo scRNAseq data (from Mohammed 2017). STAR genome built with TX1072 N-masked SNPs and "Counts_RNA_timecourse.txt" should be downloaded from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609), unzipped and stored in "./input_files/". "GENCODE_vM25_plus_Xert.gtf" should be downloaded from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".
- "3_master_analyze_RNAseq.sh": Performs analysis of RNA-seq time course and plots Fig. 3c, d-i and Extended Fig. E6b-e. "CPM_RNA_timecourse.txt" and "DESeq2_total" should be downloaded from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609) and stored in "./input_files/". "TPM_rnaseq_pacini.txt" should be retrieved from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".


All scripts that are used by the master script are stored in "./scripts/". All files that are necessary to run the master scripts that are not detailed above are stored in "./input_files/".
