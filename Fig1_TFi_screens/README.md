# Fig. 1: A CRISPR screen targeting expressed TFs detects novel Xist regulators

## Description
This folder contains all code necessary to reproduce the analysis for Figure 1 and Extended Figures E1-E3. Master scripts have to be started from the "./master" directory. The file structure of this github and files/directories retrieved from zenodo should be kept intact. All master scripts can be run independently. If not specified otherwise, any necessary files are provided in "./input_files/".


## Software dependencies and operating systems
In order to perform the analyses the following software has to be installed and available on the command line ($PATH):
- MAGeCK (v0.5.9.3) CRISPR screen analysis toolbox from "https://sourceforge.net/p/mageck/wiki/Home/"

R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- EnvStats (v2.8.1)
- flowCore (v2.10.0)
- ggcyto (v1.26.4)
- ggrepel (v0.9.2)
- ggvenn (v0.1.10)
- gridExtra (v2.3)
- opencyto (v2.10.1)
- pheatmap (v1.0.12)
- Rgb (v1.7.5)
- Rsubread (v2.12.3)
- readxl (v1.4.1)
- scales (v1.3.0)
- tidyverse (v1.3.2)
- viridis (v0.6.2)


## Reproduce analysis
All master scripts can be run independently from each other. Detailed instructions are listed below.

- "1_master_TFi_Lib_design.sh": Designs the TFi library and plots Extended Figures E1b-d. BAM files from published TT-seq data [GSE167356](github.com/EddaSchulz/Xert_paper), "TPM_rnaseq_pacini.txt" and "GENCODE_vM25_plus_Xert.gtf" should be retrieved from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424) and stored in "./input_files/".
- "2_master_TFiScreen_qc.sh": Performs read counting and guide enrichment for the TFi screen using MAGeCK. Performs quality control analysis and plots Extended Figures E1f-g, E2a-e. FASTQ.GZ files of the TFi screen should be retrieved from [GSE273070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273070) and stored in "./input_files/". The FASTQ.GZ files should be named according to the pattern: TFi_R[1/2]_[Selected/Unsorted/High/Low/Negative].fastq.gz. The plasmid library file should be named TFi_PlasmidLib.fastq.gz. FCS files from the TFi screen should be retrieved from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424), unzipped and stored in "./input_files/".
- "3_master_TFiMiniScreen_qc.sh": Designs the TFiMini library. Performs read counting and guide enrichment for the TFiMini screen using MAGeCK. Performs quality control analysis and plots Extended Figures E3a-g. FASTQ.GZ files of the TFiMini screen should be retrieved from [GSE273070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273070) and stored in "./input_files/". The FASTQ.GZ files should be named according to the pattern "TFiMini_R[1/2]_[Selected/Unsorted/High/Low/Negative].fastq.gz". The plasmid library file should be named "TFiMini_PlasmidLib.fastq.gz". FCS files from the TFiMini screen should be retrieved from [10.5281/zenodo.12822424](https://zenodo.org/records/12822424), unzipped and stored in "./input_files/".
- "4_master_TFi_screens_analysis.sh": Performs analysis of the TFi and TFiMini screens. Plots Figures 1b-d and Extended Figures E2f-h, E3h-i.

All scripts that are used by the master scripts are stored in "/./scripts/". All files that are necessary to run the master scripts that are not detailed above are stored in "./input_files/".
