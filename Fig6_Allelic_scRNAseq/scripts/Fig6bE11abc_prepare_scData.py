#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import gtfparse
import sys

# Set the working directory
wd = sys.argv[1]
input_dir = os.path.join(wd, 'input_files/')
output_dir = os.path.join(wd, 'output_files/')
os.chdir(wd)

# Get list of chrX genes
gtf_filename = os.path.join(input_dir, "GENCODE_vM25_plus_Xert_Linx.gtf")
gtf = pd.DataFrame(gtfparse.read_gtf(gtf_filename))
chrX_genes = gtf[gtf[0] == 'chrX'][8].unique() # col 0 is chromosome and col 8 is gene id



######################################################################
# Read-in data
######################################################################

### Read-in count tables
counts_all = pd.read_table(os.path.join(output_dir, 'scRNAseq_combi_counts.tsv'))
counts_ha1 = pd.read_table(os.path.join(output_dir, 'scRNAseq_combi_counts_haplo1_chrX.tsv'))
counts_ha2 = pd.read_table(os.path.join(output_dir, 'scRNAseq_combi_counts_haplo2_chrX.tsv'))

### Cell barcodes files
barcode_files_rep1 = {
	'barcodes_WT_rep1.txt': 'WT',
	'barcodes_dFtxXert_rep1.txt': 'dFtx-Xert'
	}
	
barcode_files_rep2 = {
	'barcodes_WT_rep2.txt': 'WT',
	'barcodes_dFtxXert_rep2.txt': 'dFtx-Xert'
	}

# Ensembl IDs of genes of interest to be plotted
genes_of_interest = {
    'Xist': 'ENSMUSG00000086503.4'
    }

# Define what samples will be plotted and in what order
order = ['WT','dFtx-Xert']






######################################################################
# Bring count tables into the correct format and assign samples
######################################################################

### Define function to filter count tables by cell barcodes and assign the correct samples
# Note: counts_ha1_total can either be the counts for haplotype1 or the total (non-allelic) counts
# Note: for non-allelic counts, only counts_ha1_total should be provided
# Note: for allelic counts, counts_ha1_total and counts_ha2 should be provided
# Note: for allelic counts, the two counts tables will be merged
def prepareCounts(barcode_files, replicate, counts_ha1_total, counts_ha2=None):
    
    # Initialize an empty dataframe
    counts_df = pd.DataFrame()
    
    # Loop through barcode files (i.e. the different samples)
    for barcode_file, name in barcode_files.items():
        # Read-in cell barcodes
        cell_barcodes = pd.read_table(os.path.join(input_dir, barcode_file), header=None)
        # Filter count tables based on barcodes
        counts_filtered1 = counts_ha1_total[counts_ha1_total['cell'].isin(cell_barcodes[0])] # barcodes are in column 0
        if counts_ha2 is not None: # for allelic counts
            counts_filtered2 = counts_ha2[counts_ha2['cell'].isin(cell_barcodes[0])] # barcodes are in column 0
            # Merge counts from haplotype 1 and 2 into one dataframe
            counts_merged = counts_filtered1.merge(counts_filtered2, on=['cell','gene'], how='outer', suffixes=('_ha1', '_ha2'))
        else: # for non-allelic counts
            counts_merged = counts_filtered1.copy()
            
        # Add an extra column containing sample name
        counts_merged['sample'] = name
        # Concatenate dataframes for each sample
        if counts_df.empty == True:
            counts_df = counts_merged
        else:
            counts_df = pd.concat([counts_df.reset_index(drop=True), counts_merged], axis=0, ignore_index=True)

    # Replace NAs with 0s
    counts_df = counts_df.fillna(0)
    
    # Convert values to integers
    if 'count' in counts_df.columns: # for non-allelic counts
        counts_df['count'] = counts_df['count'].astype(int)
    else: # for allelic counts
        counts_df['count_ha1'] = counts_df['count_ha1'].astype(int)
        counts_df['count_ha2'] = counts_df['count_ha2'].astype(int)
        
    # Add replicate column
    counts_df['rep'] = replicate
    
    return counts_df

### Allelic counts
# Apply function to filter count tables by cell barcodes and assign the correct samples
allele_counts_rep1 = prepareCounts(barcode_files_rep1, 'rep1', counts_ha1, counts_ha2)
allele_counts_rep2 = prepareCounts(barcode_files_rep2, 'rep2', counts_ha1, counts_ha2)

### Total (non-allelic) counts
# Apply function to filter count tables by cell barcodes and assign the correct samples
total_counts_rep1 = prepareCounts(barcode_files_rep1, 'rep1', counts_all)
total_counts_rep2 = prepareCounts(barcode_files_rep2, 'rep2', counts_all)
# Wide format of count table
total_counts_rep1_wide = total_counts_rep1.pivot_table(index=['cell', 'sample'], columns='gene', values='count')
total_counts_rep2_wide = total_counts_rep2.pivot_table(index=['cell', 'sample'], columns='gene', values='count')
# Add row indeces
total_counts_rep1_wide.reset_index(inplace=True)
total_counts_rep2_wide.reset_index(inplace=True)
# Replace NaN with 0
total_counts_rep1_wide = total_counts_rep1_wide.fillna(0)
total_counts_rep2_wide = total_counts_rep2_wide.fillna(0)







######################################################################
# Create anndata object of non-allelic counts, basic filtering and QC
######################################################################

### Define function to create an anndata object of total counts:
    # performs basic QC
    # normalizes gene expression
    # calculates chrX allelic ratios and adds info as metadata
    # adds extra metadata needed for plotting later
    # functions also returns barcodes of cells which get filtered out
# Note: takes the wide format dataframe of the total counts as input (as well as the allelic counts to calculate allelic ratios)
# Note: returns anndata object of the normalized gene expression values plus info on allelic ratio etc as metadata
def anndata_norm_totalCounts(counts_wide, allele_counts, rep):
      
    ### Initialize anndata object (using non-allelic counts)   
    # Extract the gene, sample and cell info from the count table
    gene_list = counts_wide.columns.tolist()[2:]  # exclude the first two entries because they contain 'cell' and 'sample'
    sample_list = counts_wide['sample'].tolist()
    cell_list = counts_wide['cell'].tolist()
    # Extract the count values as a matrix
    X = counts_wide.iloc[:,2:].values # exclude the first two entries because they contain 'cell' and 'sample'
    # Convert counts matrix to anndata object
    adata = ad.AnnData(X)
    # Add observations and variables to anndata object
    adata.obs_names = cell_list
    adata.var_names = gene_list
    # Add the sample info to the anndata object
    adata.obs["sample"] = sample_list
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
        
    ### Calculate allelic ratios (using allelic counts) and add as metadata
    # Remove Xist before calculating allelic ratios
    allele_counts_noXist = allele_counts[~allele_counts['gene'].str.contains(genes_of_interest['Xist'])]
    # Collapse counts per cell
    allele_counts_noXist_perCell = allele_counts_noXist.groupby(['cell', 'sample']).sum().reset_index()
    # Select only the necessary columns
    allele_counts_noXist_perCell = allele_counts_noXist_perCell[['cell', 'sample', 'count_ha1', 'count_ha2']]
    # Calculate allelic ratios
    allele_counts_noXist_perCell['allelic_ratio'] = allele_counts_noXist_perCell['count_ha1'] / (allele_counts_noXist_perCell['count_ha1'] + allele_counts_noXist_perCell['count_ha2'])
    # Add allelic ratio infos as metadata to anndata object
    adata.obs['allelic_ratio'] = allele_counts_noXist_perCell['allelic_ratio'].values
      
    ### Add extra metadata for plotting later    
    # Add the differentiation time point and the replicate to the metadata
    adata.obs['day'] = pd.Categorical(allele_counts_noXist_perCell['sample'].apply(lambda x: 'day 0' if x == '2i' else ('day 2' if x == 'NTC-0' else 'day 4')))
    adata.obs['replicate'] = rep    
    # Add raw Xist counts and label on Xist_pos or Xist_neg to metadata
    adata.obs['Xist_raw'] = pd.DataFrame(adata[:, adata.var.index == genes_of_interest['Xist']].X.flatten()).values
    adata.obs['Xist_pos_neg'] = np.where(adata.obs.Xist_raw > 0, 'Xist_pos', 'Xist_neg')
 
    ### QC  
    # Extract the ensembl IDs for mitochondrial genes
    mito_gene_names = sc.queries.mitochondrial_genes("mmusculus", attrname="ensembl_gene_id")
    adata.var["mt"] = adata.var_names.str.startswith(tuple(mito_gene_names["ensembl_gene_id"]))
    # Calculate common quality control (QC) metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    # General QC plots
    sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True)
    sc.pl.violin(adata, "pct_counts_mt", groupby="sample", jitter=0.4, multi_panel=True)
    # Filter out cells out with more than 3% mitochondrial reads
    all_cells = adata.obs_names # Save the barcodes of all cells before filtering
    adata = adata[adata.obs.pct_counts_mt < 3, :]
    # Get the barcodes of the cells that were filtered out
    filtered_cells = np.setdiff1d(all_cells, adata.obs_names)
      
    ### Normalize counts  
    # Counts are normalized to 10,000 reads per cell and log-transformed
    sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
        
    ### Add extra metadata for plotting later
    # Add the normalized expression values of the genes_of_interest (Xist) to the metadata
    for gene in genes_of_interest.keys():
        adata.obs[gene] = pd.DataFrame(adata[:, adata.var.index == genes_of_interest[gene]].X.flatten()).values

    return adata, filtered_cells

# Apply function to obtain an anndata object with normalized gene expression values
adata_rep1, barcodes_filteredCells_rep1 = anndata_norm_totalCounts(total_counts_rep1_wide, allele_counts_rep1, 'rep1')
adata_rep2, barcodes_filteredCells_rep2 = anndata_norm_totalCounts(total_counts_rep2_wide, allele_counts_rep2, 'rep2')






######################################################################
# Filter out XO cells and combine both replicates; filter for Xist+ cells
######################################################################

### Filter out XO cells from anndata object
adata_rep1_noXO = adata_rep1[~((adata_rep1.obs.Xist < 0.5) & ((adata_rep1.obs.allelic_ratio < 0.15) | (adata_rep1.obs.allelic_ratio > 0.85)))]
adata_rep2_noXO = adata_rep2[~((adata_rep2.obs.Xist < 0.5) & ((adata_rep2.obs.allelic_ratio < 0.15) | (adata_rep2.obs.allelic_ratio > 0.85)))]
# Get the barcodes of the XO cells which get filtered out
barcodes_filteredCells_XO_rep1 = np.array(adata_rep1[(adata_rep1.obs.Xist < 0.5) & ((adata_rep1.obs.allelic_ratio < 0.15) | (adata_rep1.obs.allelic_ratio > 0.85))].obs_names)
barcodes_filteredCells_XO_rep2 = np.array(adata_rep2[(adata_rep2.obs.Xist < 0.5) & ((adata_rep2.obs.allelic_ratio < 0.15) | (adata_rep2.obs.allelic_ratio > 0.85))].obs_names)

### Combine the two replicates
adata_combined = ad.concat([adata_rep1_noXO, adata_rep2_noXO], axis=0, join='outer')
# Replace NaN values with 0 in the .X matrix
adata_combined.X = np.nan_to_num(adata_combined.X)


### Filter for Xist+ cells
adata_rep1_noXO_XistPos = adata_rep1_noXO[adata_rep1_noXO.obs.Xist_pos_neg == 'Xist_pos']
adata_rep2_noXO_XistPos = adata_rep2_noXO[adata_rep2_noXO.obs.Xist_pos_neg == 'Xist_pos']
# Get the barcodes of the Xist- cells that get filtered out
barcodes_filteredCells_XistNeg_rep1 = np.array(adata_rep1_noXO[adata_rep1_noXO.obs.Xist_pos_neg == 'Xist_neg'].obs_names)
barcodes_filteredCells_XistNeg_rep2 = np.array(adata_rep2_noXO[adata_rep2_noXO.obs.Xist_pos_neg == 'Xist_neg'].obs_names)

### Combine the two replicates
adata_combined_XistPos = ad.concat([adata_rep1_noXO_XistPos, adata_rep2_noXO_XistPos], axis=0, join='outer')
# Replace NaN values with 0 in the .X matrix
adata_combined_XistPos.X = np.nan_to_num(adata_combined_XistPos.X)





######################################################################
# Plot PCA
######################################################################

### Plot PCA using total counts
def plot_pca(adata, color, palette, colormap, filename=None):
    
    # Feature selection
    sc.pp.highly_variable_genes(adata, n_top_genes=500, batch_key="sample")
    #sc.pl.highly_variable_genes(adata)
    
    # PCA
    sc.tl.pca(adata)
    #sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    fig, ax = plt.subplots()
    sc.pl.pca(adata,color=color,dimensions=(0, 1), size=40, ax=ax, palette=palette, color_map=colormap)
    plt.show()
    if filename == True:
        fig.savefig(os.path.join(output_dir, color + '_PCA.pdf'), bbox_inches='tight')


# Apply function to create PCA plots
# Note: adata_combined contains both replicates and Xist-pos and Xist-neg cells (cells have been filtered for basic QC and XO)
plot_pca(adata_combined, 'Fig6b_sample', None, None, filename=True)
plot_pca(adata_combined, 'FigE11a_replicate', None, None, filename=True)


    





######################################################################
# Allelic counts: plot number of reads and genes detected (chrX only)
######################################################################

# Filter cells from allelic count table
# Note: based on total counts cells with high mitochondrial reads and XO cells were filtered out --> filter out same cells
# Concatenate barcodes
total_filteredCells_rep1 = np.concatenate((barcodes_filteredCells_rep1, barcodes_filteredCells_XO_rep1))
total_filteredCells_rep2 = np.concatenate((barcodes_filteredCells_rep2, barcodes_filteredCells_XO_rep2))
# Filter out cells
allele_counts_rep1_filtered = allele_counts_rep1[~allele_counts_rep1['cell'].isin(total_filteredCells_rep1)]
allele_counts_rep2_filtered = allele_counts_rep2[~allele_counts_rep2['cell'].isin(total_filteredCells_rep2)]
# Concatenate replicates
allele_counts_concat = pd.concat([allele_counts_rep1_filtered, allele_counts_rep2_filtered], axis=0)
# Add the counts from haplotype 1 and haplotype 2
allele_counts_concat['count_ha1_ha2'] = allele_counts_concat['count_ha1'] + allele_counts_concat['count_ha2']
# Wide format of count table
allele_counts_concat_wide = allele_counts_concat.pivot_table(index=['cell', 'sample'], columns='gene', values='count_ha1_ha2')
# Add row indeces
allele_counts_concat_wide.reset_index(inplace=True)
# Replace NaN with 0
allele_counts_concat_wide = allele_counts_concat_wide.fillna(0)


### Define function to create an anndata object of allelic counts:
    # plots basic QC
    # outputs the number of detected genes (only chrX possible, allelic) and total allelic counts per cell
# Note: takes the wide format dataframe of the allelic counts as input (merged replicates)
def anndata_allelicCounts(allele_counts_wide, filename=None):
       
    ### Initialize anndata object (using allelic counts)
    # Extract the gene, sample and cell info from the count table
    gene_list = allele_counts_wide.columns.tolist()[2:]  # exclude the first two entries because they contain 'cell' and 'sample'
    sample_list = allele_counts_wide['sample'].tolist()
    cell_list = allele_counts_wide['cell'].tolist()
    # Extract the count values as a matrix
    X = allele_counts_wide.iloc[:,2:].values # exclude the first two entries because they contain 'cell' and 'sample'
    # Convert counts matrix to anndata object
    adata = ad.AnnData(X)
    # Add observations and variables to anndata object
    adata.obs_names = cell_list
    adata.var_names = gene_list
    # Add the sample info to the anndata object
    adata.obs["sample"] = sample_list
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
       
    ### QC   
    # Calculate common quality control (QC) metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    # General QC plots
    fig = sc.pl.violin(adata,["n_genes_by_counts", "total_counts"],jitter=0.4,multi_panel=True, show=False)
    # Save plot
    if filename == True:
        fig.savefig(os.path.join(output_dir, 'FigE11bc_NumGenesChrX_mutant_AllelicCounts_perCell.pdf'), bbox_inches='tight')
    # Calculate median values
    median_n_genes_by_counts = np.median(adata.obs['n_genes_by_counts'])
    median_total_counts = np.median(adata.obs['total_counts'])
        
    # Print median values
    print("Median n_genes_by_counts:", median_n_genes_by_counts)
    print("Median total_counts:", median_total_counts)
    

# Apply function to plot the number of chrX detected genes and allelic counts
anndata_allelicCounts(allele_counts_concat_wide, filename=True)







######################################################################
# Save tables
######################################################################

### Save normalized counts plus metadata (total counts, normalized, QC and XO filtered, Xist pos+neg)
df = adata_combined.to_df().join(adata_combined.obs)
df.to_csv(os.path.join(output_dir, 'CountTable_plusMetadata_totalCounts_mutants.txt'), sep='\t', index=True, na_rep='NaN')

# Using the allelic counts which have been filtered for high mitochondrial reads and XO cells, filter out Xist- cells
allele_counts_rep1_filtered_XistPos = allele_counts_rep1_filtered[~allele_counts_rep1_filtered['cell'].isin(barcodes_filteredCells_XistNeg_rep1)]
allele_counts_rep2_filtered_XistPos = allele_counts_rep2_filtered[~allele_counts_rep2_filtered['cell'].isin(barcodes_filteredCells_XistNeg_rep2)]
# Concatenate replicates
allele_counts_XistPos_concat = pd.concat([allele_counts_rep1_filtered_XistPos, allele_counts_rep2_filtered_XistPos], axis=0)
# Add normalized Xist expression and raw Xist counts (from total counts)
Xist_norm_raw_df = pd.DataFrame(adata_combined_XistPos.obs.Xist_raw)
Xist_norm_raw_df.columns = ["total_Xist_counts"]
Xist_norm_raw_df['Xist_norm'] = pd.DataFrame(adata_combined_XistPos.obs.Xist)['Xist']
Xist_norm_raw_df = Xist_norm_raw_df.reset_index().rename(columns={'index': 'cell'})
allele_counts_XistPos_concat_plusXist = pd.merge(allele_counts_XistPos_concat, Xist_norm_raw_df, on='cell', how='left')
# Save the table
allele_counts_XistPos_concat_plusXist.to_csv(os.path.join(output_dir, 'AllelicCounts_mutant_plusTotalXist.txt'), sep='\t', index=True, na_rep='NaN')

