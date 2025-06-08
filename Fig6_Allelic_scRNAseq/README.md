# Fig. 6: Distal REs are required for efficient silencing during random XCI

## Description
This folder contains all code necessary to reproduce the analysis for Figure 6 and Extended Figure E10. Please note that the WT/dFtx-Xert samples were sequenced together with an unrelated experiments using the MULTI-seq protocol. Provided FASTQ.GZ files (GSE273072) were filtered for MULTI-seq barcodes belonging to the experiments included in Schwämmle et al., 2025 using UMI-tools (v1.1.4) with [umi_tools extract --extract-method=string --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --whitelist=sample_whitelist.tsv --stdin total_R1.fastq.gz --read2-in=total_R2.fastq.gz --stdout subset_R1.fastq.gz --read2-out=subset_R2.fastq.gz]. Example code for the demultiplexing is provided in "./example_demultiplex/". Code to replicate the alignment from the filtered FASTQ.GZ files is instead provided in "./master/" and "./scripts/". For the XistTitration experiment scripts starting from processed count tables are provided. The pre-processing was performed analogous to the code specified for the WT/dFtx-Xert samples.  Master scripts have to be started from the "./master/" directory. The file structure of this github and files/directories retrieved from zenodo should be kept intact. All master scripts can be run independently. If not specified otherwise, any necessary files are provided in "./input_files/".


## Software dependencies and operating systems
In order to perform the analyses the following software has to be installed and available on the command line ($PATH):
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- STAR (v2.5.4a) aligner from "https://github.com/alexdobin/STAR"
- Subread (v2.0.6) suite for processing NGS read data from "https://subread.sourceforge.net/"
- UMI-tools (v1.1.4) collection of tools to deal with unique molecular identifiers and cell barcodes from "https://github.com/CGATOxford/UMI-tools"

Python (v3.10.10) script can be run with the following libraries:
- anndata (v0.10.3)
- gtfparse (v2.5.0)
- matplotlib (v3.8.2)
- numpy (v1.26.2)
- pandas (v2.1.4)
- scanpy (v1.9.4)


R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- drc (v3.0.1)
- egg (v0.4.5)
- gridExtra (v2.3)
- pheatmap (v1.0.12)
- Rgb (v1.7.5)
- tidyverse (v1.3.2)
- viridis (v0.6.2)


## Reproduce analysis
All master scripts can be run independently from each other. Detailed instructions are listed below.

- "1_master_analyze_CRISPRi_cont.sh: Processes and analyzes qPCR and FISH data from sgFtx-Xert-REs CRISPRi line. Plots Extended Fig. E10b-d. FISH images (.czi files) should be retrieved from [10.5281/zenodo.15617118](https://zenodo.org/records/15617118) and stored in "./input_files/CRISPRi_cont_FISH/".
- "2_master_process_scRNAseq.sh": Processes scRNA-seq data. Generates count tables (allelic and complete). Plots Fig. 6b and Extended Fig. E10e Gene expression FASTQ.GZ files should be retrieved from [GSE273072](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273072), renamed by the pattern "scRNAseq_rep[1|2]_[WT|dFtxXert]_R[1|2].fastq.gz" and stored in "./input_files/scRNAseq_fastq/". "TX1072_SNPs.vcf", "mm10_edit.fa", "mm10_edit.fa.fai", "GENCODE_edit.gtf" and "cellbarcode_whitelist.txt" should be retrieved from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609) and stored in "./input_files/". "GENCODE_vM25_plus_Xert_Linx.gtf" should be retrived from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".
- "3_master_analyze_scRNAseq.sh": Analyzes X-linked allelic expression in the scRNA-seq data. Plots Fig. 6c-e and Extended Fig. E10f-i. "CPM_RNA_timecourse.txt" and "cellbarcode_whitelist.txt" should be downloaded from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609) and stored in "./input_files/". "GENCODE_vM25_plus_Xert_Linx.gtf" should be downloaded from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".
- "4_master_analyze_XistTitration.sh": Analyzes Xist titration scRNA-seq data. Plots Fig. 6g-i, k-o and Extended Fig. E10j-u. "XistTitration_R1_counts.tsv" and "XistTitration_R2_counts.tsv" should be downloaded from [10.5281/zenodo.15617118](https://zenodo.org/records/15617118) and stored in "./input_files/". "CPM_RNA_timecourse.txt" should be downloaded from [10.5281/zenodo.13254609](https://zenodo.org/records/13254609) and stored in "./input_files/". "GENCODE_vM25_plus_Xert_Linx.gtf" should be downloaded from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".


All scripts that are used by the master scripts are stored in "/./scripts/". All files that are necessary to run the master scripts that are not detailed above are stored in "./input_files/".
