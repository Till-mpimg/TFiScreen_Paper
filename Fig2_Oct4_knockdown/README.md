# Fig. 2: OCT4 acts as an Xist activator during differentiation

## Description
This folder contains all code necessary to reproduce the analysis for Figure 2 and Extended Figures E5. The master script has to be started from the "./master" directory. The file structure of this github should be kept intact. Any necessary files to run the master script are provided in "./input_files/". Microscope image from RNA-FISH experiments are provided at [10.5281/zenodo.12821363](https://zenodo.org/records/12821363) (not required to reproduce the analysis detailed here).


## Software dependencies and operating systems
R scripts can be run using R (v4.2) software ("https://cran.r-project.org/"). The following R libraries are required:
- egg (v0.4.5)
- EnvStats (v2.8.1)
- gg4hx (v0.2.8)
- gridExtra (v2.3)
- readxl (v1.4.1)
- scales (v1.3.0)
- tidyverse (v1.3.2)


## Reproduce analysis
Detailed instructions are listed below.

- "1_master_analyze_Oct4_knockdown.sh": Analyzes qPCR data from differentiation experiments. Analyzes qPCR and RNA-FISH data from Oct4 knockdown experiments. Plots Figures 2b-c and Extended Figures E4a-d.
- BIGWIG tracks and NARROWPEAK files of OCT4 binding in male ESCs and EpiLCs were generated as described in [github.com/EddaSchulz/XertPaper](https://github.com/EddaSchulz/Xert_paper/tree/main/NGS_alignment) with master script "(x) master_Buecker_Stadler_ChIPseq.sh" and can be retrieved from [10.5281/zenodo.15617118](https://zenodo.org/records/15617118).

All scripts that are used by the master script are stored in "./scripts/". All files that are necessary to run the master script are stored in "./input_files/".
