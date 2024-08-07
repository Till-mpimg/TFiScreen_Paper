#####################################################################/
# ~~~~ sce multi analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Jonathan Froehlich
# Schulz lab, Berlin
#
# edited by Till Schwämmle
#####################################################################/
# This script serves as example code for the filtering and MULTI barcode association for the single cell experiment published in Schwämmle et al. 2024
# Filtered FASTQ.GZ files with the WT/dFtx-Xert samples are provided at GSE273072
# The libraries were prepared together with an unrelated experiment
# This code exemplifies how the files were split and how initial quality filtering was performed
# Split lists of the final barcodes per sample are provided in "./Fig6_Allelic_scRNAseq/input_files" and used for the master scripts. 

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
## Manual outline ####
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# 1. Setup
# Define paths
# Load libraries
# Define ggplot theme

# 2. Load, normalize, QC, filter
# Start loop to process each of the sc datasets (e.g. JJF002, JJF003, JJF004)
#     Load count matrix, reformat sce (SingleCellExperiment object) slightly
#     Process MULTI-seq fastqs with deMULTIplex2
#     Add MULTI-seq annotations to cells
#     Normalize by library size (currently simple size factors from library size)
#     log2 transform
#     Feature selection (currently top 10% highly variable genes)
#     PCA, Louvain clustering, UMAP
#     Add QC metrics (mito, etc.)
#
#     "sce_overview_plots" function:
#     - number of cells per sample. 
#     - text summarizing some metrics (e.g. cells filtered by QC, median UMIs per cell, etc.)
#     - PCA plots (e.g. by cluster label, by sample name)
#     - UMAP plots (e.g. by cluster label, by sample name)
#     - violin plots (e.g. UMI counts, detected genes, mito%)
#
#     unfiltered sce: 
#         - Output: pdf from sce_overview_plots
#         - Output: export sce object as RDS file
# 
#     QC filtered sce (e.g. mito):
#         Normalization and transformation again (because size factor normalization can be affected by low-quality cells)
#         Feature selection again
#         PCA, Louvain clustering, UMAP again
#         - Output: pdf from sce_overview_plots
#         - Output: export sce object as RDS file
#
#     Singlet filtered sce:
#         Opt. Filter for titration samples only (remove Tills samples)
#         Normalization and transformation again (because size factor normalization can be affected by low-quality cells)
#         Feature selection again
#         PCA, Louvain clustering, UMAP again
#         - Output: pdf from sce_overview_plots
#         - Output: export sce object as RDS file
#         - Output: export list of cell barcodes as csv
#         - Output: export the cell barcodes with MULTI-barcode annotation (e.g. sample assignment)
#
#     => Output of the loop, for each dataset: 
#     - 3 sce objects (unfiltered, after QC filter, after singlet filter)
#     - 3 pdf files with overview plots (unfiltered, after QC filter, after singlet filter)
#     - table with cell barcodes (csv)
#     - table with cell barcodes and MULTI-seq annotations (csv)

# 3. DGE
# Read and combine sce files from each dataset (JJF002, JJF003, JJF004)
# Aggregate gene counts by treatment
# Convert to data.frame, wrangle slightly
# Run DESeq
# List of genes to highlight, using several filters (e.g. all genes with padj < 0.05, top and bottom 25 genes by fold change, minimum counts)
# - Create volcano plots with highlighted genes, others with metadata overlaid (e.g. Oct4 targets, Oct4 motif, etc.)

# 4. Combine sce objects, just to make plotting easier later
# Read and combine sce files from each dataset (JJF002, JJF003, JJF004), "combined_sce"

# 5. Overview plots of dataset metrics (e.g. UMIs/cell)
# - Unfiltered sce: Cell Counts by Droplet Type and Dataset  TODO solve error
# - Filtered sce, my samples: Cell Counts distribution
# - Filtered sce, my samples: UMIs per Cell, Genes per Cell
# - Filtered sce, for each dataset: UMIs per Cell, Genes per Cell

# 6. Plot Genes expression
# Here we plot the expression of genes of interest, e.g. Oct4, Xist, transgene, and the genes from DGE analysis
# - Violin plots, for each assay type ("counts", "logcounts", "shifted_log"), separately for each dataset
# - Violin plots, for each assay type ("counts", "logcounts", "shifted_log"), all datasets next to each other (=grouped by dataset_id)
# - Scatter plots ("counts", "logcounts") for some pairs of interesting genes (e.g. Oct4 x transgene, Xist x Tsix)

# 7. Binning by Oct4 expression
# Binning by Oct4 expression, for each assay type ("counts", "logcounts", "shifted_log"), by size and by width, separately for each dataset
# Take each sample and bin by Oct4 expression, not so useful because within every sample there are cells that express Oct4 very high or low

# 8. Plot means from three datasets
# Here we plot the mean expression of genes of interest, e.g. Oct4, Xist, transgene, and the genes from DGE analysis
# Either by sample_name, or by the differently calculated Oct4 bins
# For each assay type ("counts", "logcounts", "shifted_log")
# One .pdf file for each binning approach and assay type, and each assay type of gene expression (different combinations)
# - Dot plots of mean values from 3 datasets plus their overall mean

# 9. On unfiltered_sce: binning by Oct4 expression
# ...
# 10. On unfiltered_sce: plot means from three datasets
# ...

# 11. Heatmap
# Heatmap of gene expression
# in progress

# Retired
# 12. Plot library size (=total UMIs/cell) across several sce objects saved as .rds, to see if library size is different between samples (which have different cell volumes)
# 13. Plot median UMIs per cell across several sce objects (e.g. different reference genomes), to see what effect the reference genome has on the median UMIs per cell
# 14. in progress: Histograms of gene expression levels (to visualize normalization steps), was too slow at some point
# 15. Plot gene expression, earlier versions of the plots from script v0.1.0
# 16. Find markers (alternative to DESeq, working on one replicate)


###################################################################################################################/
# Define paths ####
###################################################################################################################/

# set working directory: 
working_directory <- "PATH/TO/WORK_DIRECTORY"

# define sub directories (will be created if dont exist) 
data_dir <- "data/"
output_dir <- "output/"

# settings table to run script over several datasets 
# contains file paths to fastq files
path_settings_table <- "/project/agsgpa/Jonathan/datasets/scRNAseq/20231206_JJF002-004_Oct4/data/settings_20240202_JJF002-004_Oct4.xlsx"
sheet_settings_table <- "file_paths"


###################################################################################################################/
# Define paths If running code on only one dataset ####
###################################################################################################################/

# # file paths and settings
# # count matrix output of STARsolo
path_sc_data <- "PATH/TO/starsolo_diploidGenome_mm10-scLibrary/Solo.out/GeneFull_Ex50pAS/filtered/"
run_id <- "RUN_ID"

# # info for MULTI-seq demultiplex2 package
# # MULTI-seq fastq.gz files
prefix_multi_fastq_names <- "Library-MULTIseq"
path_multi_data1 <- "PATH/TO/Library-MULTIseq-R1.fastq.gz"
path_multi_data2 <- "PATH/TO/Library-MULTIseq-R2.fastq.gz"
# # table which includes which sample was barcoded with which MULTI-seq barcode
path_multi_to_sample_map <- "FOUND/AT/Fig6_Allelic_scRNAseq/example_demultplex/MULTI-seq_barcode_to_sample_map.xlsx"
sheet_multi_to_sample_map <- "20231122-seqLibrary" |  "20231124-seqLibrary" #specify sheet in excel table
# # good thresholds 
demultiplex2_cell_threshold <- 200 #ADAPT TO LIBRARY
demultiplex2_tag_threshold <- 2000 #ADAPT TO LIBRARY
# # MULTI lib total fragments
MULTI_lib_fragments <- 67000000 #ADAPT TO LIBRARY


###################################################################################################################/
# 0 Setup ####
###################################################################################################################/

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
## Create directories ####
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# creates sub directories
if (!dir.exists(data_dir)) {dir.create(data_dir)}
if (!dir.exists(output_dir)) {dir.create(output_dir)}


###################################################################################################################/
# 1 Load libraries ####
###################################################################################################################/

# set working directory 
setwd(working_directory) 

# load libraries
# single cell
library(DropletUtils)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(scran)
library(scater)
library(uwot)
library(bluster)
library(transformGamPoi)

# deMULTIplex2
library(deMULTIplex2)
library(Cairo)

# data
library(readxl)
library(tidyverse)
library(Hmisc)

# seed for reproducibility
set.seed(1234)

# re-assign conflicted functions
select <- dplyr::select

###################################################################################################################/
# 2 Process scRNAseq counts (load, normalize, QC, filter...) ####
###################################################################################################################/

#-----------------------------------------------#
# start of the function and loop to process the sc data, performs QC filtering and make overview plots
#-----------------------------------------------#
process_sce_datasets <- function(run_id, path_sc_data, path_multi_data1, path_multi_data2, prefix_multi_fastq_names, path_multi_to_sample_map, sheet_multi_to_sample_map, demultiplex2_cell_threshold, demultiplex2_tag_threshold) {

  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Load count matrix ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  
  # Here we load the "filtered" count matrix from STARsolo. This is already filtered for empty droplets. 
  # Alternatively one could use the "raw" output and filter in R for real cells (removing "empty" droplets) with library(DropletUtils) emptyDrops, see https://bioconductor.org/books/3.15/OSCA.advanced/droplet-processing.html#qc-droplets
  sce <- read10xCounts(path_sc_data)

  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Reformat sce slightly  ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  # Fix duplicate gene symbols
  # Check for duplicate symbols (I have encountered a few, e.g. 50, once, might have been specific to the "optimized" reference genome of Pool et al.)
  duplicate_symbols <- duplicated(rowData(sce)$Symbol) | duplicated(rowData(sce)$Symbol, fromLast = TRUE)
  # Make the symbols unique (e.g. adding a ".1" to the end of the symbol)
  rowData(sce)$Symbol[duplicate_symbols] <- make.unique(rowData(sce)$Symbol[duplicate_symbols])
  # Assign unique gene symbols as row names
  rownames(sce) <- rowData(sce)$Symbol
  
  # Use gene symbol instead of gene ID
  rowData(sce)$Symbol <- as.character(rowData(sce)$Symbol) # Convert Symbol to character if it's not already
  rownames(sce) <- rowData(sce)$Symbol
  
  # Add cell_barcode to colData
  sce$cell_barcode <- gsub("-1", "", sce$Barcode) 
  
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Process MULTI-seq fastqs with deMULTIplex2  ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # https://github.com/Gartner-Lab/deMULTIplex2
  
  # Load the table that maps the MULTI-seq barcodes to the sample names
  multi_to_sample_df <- as.data.frame(read_xlsx(path_multi_to_sample_map, sheet_multi_to_sample_map))
  
  # Modify the key column, to remove underscores (e.g. "A_12" -> "A12")
  multi_to_sample_df$well_rosen_2022_modified <- gsub("_", "", multi_to_sample_df$well_rosen_2022)

  # Once the below section was run, it can be skipped and the assignments can be loaded from file
  ### ----------------------- MARKER: comment from here in case of skipping the MULTI-seq demultiplex2 step ####

  # Copy fastq.gz files from server to local directory (otherwise error because too many other fastq in same folder)
  file.copy(path_multi_data1, data_dir)
  file.copy(path_multi_data2, data_dir)
  
  read_table <- readTags(dir = data_dir,
                         name = prefix_multi_fastq_names, # prefix of fastq file names
                         barcode.type = "MULTIseq",
                         assay = "RNA")
  
  # Delete the files again
  file.remove(file.path(data_dir, basename(path_multi_data1)))
  file.remove(file.path(data_dir, basename(path_multi_data2)))

  # Current MULTI-seq BC oligo sequences, provided with the package
  data(multiseq_oligos) 
  tag_ref <- deMULTIplex2::multiseq_oligos

  # only keep barcodes used in this experiment, listed in multi_to_sample_df$well_rosen_2022_modified
  tag_ref <- tag_ref[names(tag_ref) %in% multi_to_sample_df$well_rosen_2022_modified]
  
  # align tags
  # resulting "tag_mtx": columns = barcodes, rows = cell barcodes, values = UMI counts of barcode
  tag_mtx <- alignTags(read_table, tag_ref,  max.dist = 2)

  # filter tag_mtx for cell barcodes that are in sce (i.e. contains real cell barcodes - "empty drops" cells already removed by STARsolo/Cellranger)
  tag_mtx <- tag_mtx[rownames(tag_mtx) %in% sce$cell_barcode, ]
  
  # Core function of the deMULTIplex2 package
  # Will output list
  res <- demultiplexTags(tag_mtx,
                         plot.path = output_dir,
                         plot.name = paste0(run_id, "_"),
                         plot.diagnostics = TRUE)
  
  # export assignments
  bc_anno_df <- cbind(cell_barcode = rownames(res$assign_table), res$assign_table)
  write.csv(bc_anno_df, paste0(output_dir, paste0("/", run_id, "_MULTI_assign_table.csv")), row.names = FALSE)
  
  
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Add MULTI-seq annotations to sce  ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  
  # Add sample_name to bc_anno_df
  # Merge the dataframes
  bc_anno_df <- merge(bc_anno_df, 
                      multi_to_sample_df[, c("well_rosen_2022_modified", "sample_name")], 
                      by.x = "final_assign", 
                      by.y = "well_rosen_2022_modified",
                      all.x = TRUE)

  # Replace NA in sample_name column with the value from final_assign
  bc_anno_df$sample_name[is.na(bc_anno_df$sample_name)] <- bc_anno_df$final_assign[is.na(bc_anno_df$sample_name)]
  
  # Merge with bc-to-sample info from MULTIplex2
  rownames(bc_anno_df) <- bc_anno_df$cell_barcode
  
  # Get the matching indices
  matched_indices <- match(colData(sce)$cell_barcode, bc_anno_df$cell_barcode)
  
  # For each column you want to add from bc_anno_df to sce, do the following:
  cols_to_add <- c("barcode_assign", 
                   "barcode_count", 
                   "droplet_type", 
                   "final_assign", 
                   "sample_name")
  
  for (col in cols_to_add) {
    # Use matched_indices to map the values from bc_anno_df to the correct rows in sce
    colData(sce)[[col]] <- bc_anno_df[matched_indices, col]
  }
  
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Normalize & Dimension reduction ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  
  # scale simple
  # Library size normalization. Library size defined as the total sum of counts across all genes for each cell.
  # The “library size factor” for each cell is then directly proportional to its library size where the proportionality constant is defined such that the mean size factor across all cells is equal to 1.
  # https://bioconductor.org/books/3.15/OSCA.basic/normalization.html#library-size-normalization
  sizeFactors(sce) <- librarySizeFactors(sce)

  # standard log2 transform 
  # some info: "Internally ... the logNormCounts() function will use log1p. This has the advantage of preserving sparsity compared to the naive log(x + 1) method, which goes through a non-sparse intermediate x + 1. For all other settings of pseudo.count, we go through the usual log(x + pseudo.count), because we're going to lose sparsity anyway. Now, log(x + 1) is a natural log, so we divide the result by log(2) to obtain log2-transformed values." If another log is desired (e.g., log10), set log=FALSE, and log-transform it afterwards manually. https://github.com/LTLA/scuttle/issues/12
  sce <- logNormCounts(sce, transform = "log", pseudo.count = 1) # transform = "log", pseudo.count = 1 are default
  
    
  # Feature selection.
  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, prop = 0.1) # "prop = 0.1": top 10% of genes with the highest biological components, HVGs ("highly variable genes")
  
  # PCA
  sce <- runPCA(sce, ncomponents = 25, subset_row = hvg) # on HVGs
  
  # Clustering
  colLabels(sce) <- clusterCells(sce, use.dimred = "PCA", BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))    
  
  # UMAP
  sce <- runUMAP(sce, dimred = "PCA")
  
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Add QC metrics ####
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  
  # Quality control using mitochondrial genes.
  
  # get mitochondria genes via gene symbol
  is_mito <- grepl("^mt-", rowData(sce)$Symbol)
  # add mito QC stats to colData
  sce <- addPerCellQCMetrics(sce, subsets = list(Mito = is_mito))
  colnames(colData(sce))

  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  ## Filter for high quality ###
  ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  
  # For counting how many cells are filtered out at each step
  # Initial number of cells
  count_initial <- ncol(sce)
  
  # Export (save) sce object as RDS file before filtering
  saveRDS(sce, paste0(output_dir, paste0("/", run_id, "_sce_unfiltered.rds")))
  
  # QC filter (e.g. mito)
  # from scuttle documentation: "This identifies low-quality cells as those that are low outliers for the log-total count, 
  # low outliers for the log-number of detected features, or high outliers for the percentage of counts in specified gene sets (e.g., mitochondrial genes, spike-in transcripts)."
  qcstats <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))
  filtered <- quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
  sum(!filtered$discard)
  sum(filtered$discard)
  sce <- sce[, !filtered$discard]

  # filter out any cells with a mitochondrial gene percentage greater than 7%
  sce <- sce[, sce$subsets_Mito_percent < 7]
  
  # For counting how many cells are filtered out at each step
  count_post_qc <- sum(!filtered$discard)  # Count after QC filter
  removed_by_qc <- count_initial - count_post_qc
  
  # Normalization and transformation
  # Done again because size factor normalization can be affected by low-quality cells (https://www.bioconductor.org/packages/release/bioc/vignettes/scuttle/inst/doc/norm.html "Low-quality cells with few expressed genes can often have negative size factor estimates").
  # Should overwrite the values from the previous normalization and transformation
  
  # scale simple
  # Library size normalization. Library size defined as the total sum of counts across all genes for each cell.
  # The “library size factor” for each cell is then directly proportional to its library size where the proportionality constant is defined such that the mean size factor across all cells is equal to 1.
  # https://bioconductor.org/books/3.15/OSCA.basic/normalization.html#library-size-normalization
  sizeFactors(sce) <- librarySizeFactors(sce)

  # standard log2 transform 
  # some info: "Internally ... the logNormCounts() function will use log1p. This has the advantage of preserving sparsity compared to the naive log(x + 1) method, which goes through a non-sparse intermediate x + 1. For all other settings of pseudo.count, we go through the usual log(x + pseudo.count), because we're going to lose sparsity anyway. Now, log(x + 1) is a natural log, so we divide the result by log(2) to obtain log2-transformed values." If another log is desired (e.g., log10), set log=FALSE, and log-transform it afterwards manually. https://github.com/LTLA/scuttle/issues/12
  sce <- logNormCounts(sce, transform = "log", pseudo.count = 1) # transform = "log", pseudo.count = 1 are default
  
  
  # Feature selection.
  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, prop = 0.1) # "prop = 0.1": top 10% of genes with the highest biological components, HVGs ("highly variable genes"), this should be around 500 genes(because we detected around 5,000 genes/cell)

  # PCA
  sce <- runPCA(sce, ncomponents = 25, subset_row = hvg) # on HVGs
  # sce <- runPCA(sce, exprs_values = "logcounts") # on all genes
  # Clustering
  # One could also increase the louvain resolution parameter to get more clusters, which might help when trying to manually assign clusters to remove till_ samples
  colLabels(sce) <- clusterCells(sce, use.dimred = "PCA", BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))    
  # UMAP
  sce <- runUMAP(sce, dimred = "PCA")
  
  # take only singlets
  sce_singlets <- sce[, colData(sce)$droplet_type == "singlet" & !is.na(colData(sce)$droplet_type)]
  sce <- sce_singlets
  
  # For counting how many cells are filtered out at each step
  count_singlets <- ncol(sce_singlets)  # Count of singlets
  removed_non_singlets <- count_post_qc - count_singlets
  
  # export list of cell barcodes as csv
  # to use for read filtering in allele-specific XCI analysis
  # Extract cell barcodes from sce
  cell_barcodes <- colData(sce)$cell_barcode
  # Write cell barcodes to a csv file
  write.csv(cell_barcodes, paste0(output_dir, run_id, "_cell_barcodes_after_filtering.csv"), row.names = FALSE)
  # Write the cell barcodes with MULTI-barcode annotation (e.g. sample assignment)
  write.csv(colData(sce), paste0(output_dir, run_id, "_cell_barcodes_with_MULTI_annotation_after_filtering.csv"), row.names = FALSE)}

 
 
 process_sce_datasets(run_id, 
                         path_sc_data, 
                         path_multi_data1, 
                         path_multi_data2, 
                         prefix_multi_fastq_names, 
                         path_multi_to_sample_map, 
                         sheet_multi_to_sample_map, 
                         as.numeric(demultiplex2_cell_threshold), 
                         as.numeric(demultiplex2_tag_threshold))
  )
