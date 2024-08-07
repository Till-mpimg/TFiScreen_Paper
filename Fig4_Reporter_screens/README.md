# Fig. 4: Reporter screens identify functional TF-RE interactions.

## Description
This folder contains all code necessary to reproduce the analysis for Figure 4 and Extended Figures E7-E9. The master script has to be started from the "./master" directory. The file structure of this github should be kept intact. Any necessary files to run the master script are provided in "./input_files/" or are detailed below.


## Software dependencies and operating systems
In order to perform the analyses the following software has to be installed and available on the command line ($PATH):
- FIMO (v5.1.1) motif enrichment analysis tool from "https://meme-suite.org/meme/"
- MAGeCK (v0.5.9.3) CRISPR screen analysis toolbox from "https://sourceforge.net/p/mageck/wiki/Home/"
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"

R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- EnvStats (v2.8.1)
- fgsea (v1.24.0)
- flowCore (v2.10.0)
- ggcyto (v1.26.4)
- gridExtra (v2.3)
- opencyto (v2.10.1)
- pheatmap (v1.0.12)
- Rsubread (v2.12.3)
- scales (v1.3.0)
- tidyverse (v1.3.2)
- viridis (v0.6.2)


## Reproduce analysis
All master scripts can be run independently from each other. Detailed instructions are listed below.

- "1_master_analyze ReporterTimecourse.sh": Analyzes FACS data from the reporter timecourse experiment. Plots Fig. 4c and Extended Fig. E7b-e. "FCS_reporter_timecourse.tar.gz" and "Gjaltema_CnT_ATAC_bam.tar.gz" should be retrieved from 10.5281/zenodo.12822424, unzipped and stored in "./input_files/".
- "2_master_process_reporter_screens.sh": Quantifies reporter screens and calculates interactions scores. Performs quality control analysis and plots Fig. 4f,g and Extended Fig. E8c-e, E9a-b. Reporter screen FASTQ.GZ files should be retrieved from GSE273069, named according to the pattern "XX_D2_REPi_[noRE|RE57L|RE57M|RE57R|RE58|RE61|RE85|RE93|RE95|RE96|RE97|RE127]_[GFP-High|GFP-Low|Unsorted]_R[1|2|3]_[1|2].fastq.gz" and stored in "./input_files/REPi_fastq/". "REPi_screen_FCS.tar.gz" should be downloaded from 10.5281/zenodo.12822424, unzipped and stored in "./input_files/".
- "3_master_analyze_interaction_scores.sh": Performs downstream analysis of reporter screens and plots Fig. 3e, h-i and Extended Fig. E9c-j. "CPM_RNA_timecourse.txt" should be downloaded from 10.5281/zenodo.13254609 and stored in "./input_files/". "mm10.fa" and "mm10.fa.fai" should be retrieved from doi.org/10.5281/zenodo.12822424 and stored in "./input_files/".

All scripts that are used by the master script are stored in "./scripts/". All files that are necessary to run the master scripts that are not detailed above are stored in "./input_files/".
