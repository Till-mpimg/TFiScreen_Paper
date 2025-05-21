# Reporter CRISPR screens decipher cis- and trans-regulatory principles at the Xist locus
Till Schwämmle, Gemma Noviello, Eleni Kanata, Jonathan J. Froehlich, Melissa Bothe, Aybuge Altay, Jade Scouarnec, Vivi-Yun Feng, Martin Vingron, Edda G. Schulz

## Abstract
Developmental genes are controlled by an ensemble of cis-acting regulatory elements (REs), which in turn respond to multiple trans-acting transcription factors (TFs). Understanding how a cis-regulatory landscape integrates information from many dynamically expressed TFs has remained a challenge. Here we apply pooled CRISPR screens to the developmentally regulated Xist locus, which is specifically activated in females to initiate X-chromosome inactivation, using both the endogenous Xist RNA and RE-reporters as screen readouts. This approach enables us to comprehensively identify Xist-controlling TFs and map their TF-RE wiring. We find a group of transiently expressed TFs that regulate proximal REs, driving the binary activation of Xist expression. These basal activators are more highly expressed in cells with two X chromosomes, potentially driving female-specific Xist upregulation. A second set of developmental TFs is upregulated later during differentiation and targets distal REs. This regulatory axis is crucial to achieve high levels of Xist RNA, which is necessary for X-chromosome inactivation. Our findings support a model for developmental gene regulation in which factors targeting proximal REs drive binary ON-OFF decisions, while factors interacting with distal REs control the transcription output.

## Description

This github repository details data and code used to perform computational analyses and data visualization. All required data is either included on the github, stored at GEO (GSE273068-GSE273072) or stored at zenodo ([10.5281/zenodo.12821363](url), [10.5281/zenodo.12821095](url), [10.5281/zenodo.12822424](url) and [10.5281/zenodo.13254609](https://zenodo.org/records/13254609)). 
The scripts are structured by figures and can be executed from master scripts that combine similar analyses (stored in "master"). All master scripts can be run independently from each other. To succesfully reproduce the analyses, the file tree structure and file names should be kept intact. The utilized versions of any packages are detailed in the individual README of the different sections. 

Fig1_TFi_screens: Contains instructions, code and files to reproduce library design and screen analysis for CRISPRi screens identifying TF regulators of Xist. Plots Fig. 1 and Extended Fig. E1-E4. (Code by Till Schwämmle)

Fig2_Oct4_knockdown: Contains instructions, code and files to reproduce analysis of Oct4 knockdown experiments. Plots Fig. 2 and Extended Fig. E5. (Code by Till Schwämmle)

Fig3_RNAseq_timecourse: Contains instructions, code and files to reproduce analysis of RNA-seq time course experiments. Plots Fig. 3 and Extended Fig. E6. (Code by Till Schwämmle and Aybuge Altay)

Fig4_Reporter_screens: Contains instructions, code and files to reproduce analysis of Reporter CRISPR screen experiments. Plots Fig. 4 and Extended Fig. E7-9. (Code by Till Schwämmle)

Fig5_Xist_activator_KDs: Contains instructions, code and files to reproduce analysis of Xist activator knockdown experiments. Plots Fig. 5 and Extended Fig. E10. (Code by Till Schwämmle)

Fig6_Allelic_scRNAseq: Contains instructions, code and files to reproduce analysis of scRNA-seq experiment, studying the effects of low Xist levels on XCI. Plots Fig. 6 and Extended Fig. E11. (Code by Jonathan J. Froehlich, Mellisa Bothe and Till Schwämmle)
