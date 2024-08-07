# Fig. 5: Basal Xist activation and full transcript levels are governed by distinct mechanisms.

## Description
This folder contains all code necessary to reproduce the analysis for Figure 5 and Extended Figures E10. The master script has to be started from the "./master" directory. The file structure of this github should be kept intact. Any necessary files to run the master script are provided in "./input_files/" or are detailed below.


## Software dependencies and operating systems
Python (v3.10.10) script can be run with the following libraries:
- czifile (v2019.7.2)
- matplotlib (v3.8.2)
- numpy (v1.26.2)
- pandas (v2.1.4)
- scikit-image (v0.22)
- scipy (v1.11.4)
- seaborn (v0.13.0) 
- sklearn (v1.3.2)

In order to perform the analyses the following software has to be installed and available on the command line ($PATH):
- Bedtools (v2.29.2) collection of C++ scripts from "https://bedtools.readthedocs.io/en/latest/"
- Bowtie2 (v2.3.5.1) aligner from "https://github.com/BenLangmead/bowtie2"
- Deeptools2 (v3.4.1) collection of PYTHON scripts from "https://deeptools.readthedocs.io/en/develop/"
- MACS2 (v2.2.7.1) PYTHON script from "https://github.com/macs3-project/MACS"
- Picard tools (v2.18.25) collection of JAVA scripts from "https://broadinstitute.github.io/picard/"
- Samtools (v1.10) collection of C scripts from "http://www.htslib.org/"
- TrimGalore (v0.6.4) PERL script from "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/"

R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- EnvStats (v2.8.1)
- gridExtra (v2.3)
- pheatmap (v1.0.12)
- readxl (v1.4.1)
- Rsubread (v2.12.3)
- scales (v1.3.0)
- tidyverse (v1.3.2)


## Reproduce analysis
All master scripts can be run independently from each other. Detailed instructions are listed below.

- "1_master_Analyze_qPCR.sh": Analyzes qPCR data following Xist activator knockdown. Plots Extended Fig. E10a-c.
- "2_master_XistAct_CutnTAG.sh": Analyzes CUT&Tag data following Xist activator knockdwon. Produces BIGWIG tracks and plots Fig. 5c. "bowtie2_mm10.tar.gz" should be retrieved from 10.5281/zenodo.13254609, unzipped and stored in "./input_files/". FASTQ.GZ files should be retrieved from GSE273068. named according to the pattern "XX_D2_H3K[27ac|4me3|4me1]_sg[NT|Zic3|Nfrkb|Otx2|Foxd3]_R[1|2|3]_[1|2].fastq.gz" and stored in "./input_files/CnT_fastq/".
- "3_master_XistAct_FISH.sh": Processes FISH images and quantifies frequency/intensity of Xist signals. Plots Fig. 5d-g and Extended Fig. E10d. FACS data should be retrieved from 10.5281/zenodo.12821095, unzipped and stored in "./input_files/czi_XistAct/".

All scripts that are used by the master script are stored in "./scripts/". All files that are necessary to run the master scripts that are not detailed above are stored in "./input_files/".
