#!/usr/bin/env python3

import os
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import gtfparse

# Set the working directory
wd = sys.argv[1]
input_dir = os.path.join(wd, 'input_files/')
output_dir = os.path.join(wd, 'output_files/')
os.chdir(wd)

# Get list of chrX genes
gtf_filename = os.path.join(input_dir, "GENCODE_vM25_plus_Xert_Linx.gtf")
gtf = pd.DataFrame(gtfparse.read_gtf(gtf_filename))
chrX_genes = gtf[gtf[0] == 'chrX'][10].unique() # col 0 is chromosome and col 8 is gene id


######################################################################
# Read-in data
######################################################################

### Read-in count tables
counts = pd.read_table(os.path.join(input_dir, 'GSE151009_UMICountMatrix.txt'))
counts_b6 = pd.read_table(os.path.join(input_dir, 'GSE151009_B6_UMICountMatrix.txt'))
counts_cast = pd.read_table(os.path.join(input_dir, 'GSE151009_Cast_UMICountMatrix.txt'))

### Reformat total counts
counts_T = counts.T
new_index = counts_T.index
new_index.columns = ['Cell']
counts_T.index = new_index

### Reformat allelic counts
counts_b6_T = counts_b6.T
counts_cast_T = counts_cast.T
new_index_b6 = counts_b6_T.index
new_index_b6.columns = ['Cell']
new_index_cast = counts_cast_T.index
new_index_cast.columns = ['Cell']
counts_b6_T.index = new_index_b6
counts_cast_T.index = new_index_cast

# Extract chrX genes for allelic counts
chrX_genes_b6 = [gene for gene in chrX_genes if gene in counts_b6_T.columns]
counts_b6_chrX = counts_b6_T[chrX_genes_b6]
chrX_genes_cast = [gene for gene in chrX_genes if gene in counts_cast_T.columns]
counts_cast_chrX = counts_cast_T[chrX_genes_cast]



######################################################################
# Create anndata object, QC, normalize, calculate allelic ratio
######################################################################

### Define function to create an anndata object of total counts:
    # performs basic QC
    # normalizes gene expression
    # calculates chrX allelic ratios and adds info as metadata
    # adds extra metadata needed for plotting later
    # functions also returns barcodes of cells which get filtered out
# Note: takes the wide format dataframe of the total counts as input (as well as the allelic counts to calculate allelic ratios)
# Note: returns anndata object of the normalized gene expression values plus info on allelic ratio etc as metadata
def anndata_norm_totalCounts(counts, counts_cast, counts_b6):
       
    ### Initialize anndata object (using non-allelic counts)  
    # Extract the gene, sample and cell info from the count table
    gene_list = counts.columns.tolist()
    cell_list = counts.index.tolist()
    # Extract the count values as a matrix
    X = counts.iloc[:].values # exclude the first two entries because they contain 'cell' and 'sample'
    # Convert counts matrix to anndata object
    adata = ad.AnnData(X)
    # Add observations and variables to anndata object
    adata.obs_names = cell_list
    adata.var_names = gene_list
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
        
    ### Calculate allelic ratios (using allelic counts) and add as metadata    
    # Remove Xist before calculating allelic ratios
    counts_cast_noXist = counts_cast.drop(columns=['Xist'])
    counts_b6_noXist = counts_b6.drop(columns=['Xist'])
    # Sum all chrX allelic counts per cell
    sum_allelic_counts_cast = counts_cast_noXist.sum(axis=1)
    sum_allelic_counts_b6 = counts_b6_noXist.sum(axis=1)
    # Concat
    sum_allelic_counts = pd.concat([sum_allelic_counts_cast, sum_allelic_counts_b6], axis=1, keys=['counts_cast', 'counts_b6'])
    # Add columns of Xist counts
    sum_allelic_counts['Xist_cast'] = counts_cast['Xist']
    sum_allelic_counts['Xist_b6'] = counts_b6['Xist'] 
    # Calculate allelic ratios
    sum_allelic_counts['allelic_ratio'] = sum_allelic_counts['counts_cast'] / (sum_allelic_counts['counts_cast'] + sum_allelic_counts['counts_b6'])    
    
    ### Add extra metadata for plotting later    
    # Add allelic ratio infos as metadata to anndata object
    adata.obs['allelic_ratio'] = sum_allelic_counts['allelic_ratio'].values    
    # Add raw Xist counts and label on Xist_pos or Xist_neg to metadata
    adata.obs['total_Xist_counts'] = counts['Xist'].values
    adata.obs['Xist_pos_neg'] = np.where(adata.obs.total_Xist_counts > 0, 'Xist_pos', 'Xist_neg')

    ### QC
    # Extract the ensembl IDs for mitochondrial genes
    mito_gene_names = sc.queries.mitochondrial_genes("mmusculus", attrname="external_gene_name")
    adata.var["mt"] = adata.var_names.str.startswith(tuple(mito_gene_names["external_gene_name"]))
    # Calculate common quality control (QC) metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    # General QC plots
    sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True)
    # Filter out cells out with more than 10% mitochondrial reads and total UMI counts <50000
    all_cells = adata.obs_names # Save the barcodes of all cells before filtering
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    adata = adata[adata.obs.total_counts > 50000, :]
    ### Filter out XO cells from anndata object
    adata = adata[~((adata.obs.total_Xist_counts < 1) & ((adata.obs.allelic_ratio < 0.2) | (adata.obs.allelic_ratio > 0.8)))]
    # Plot again
    sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True)
    
    # Get the barcodes of the cells that were filtered out
    filtered_cells = np.setdiff1d(all_cells, adata.obs_names)

    ### Normalize counts    
    # Counts are normalized to 10,000 reads per cell and log-transformed
    sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    # Add norm. Xist as metadata
    adata.obs["Xist_norm"] = pd.DataFrame(adata[:, adata.var.index == "Xist"].X.flatten()).values

    return adata, filtered_cells

adata, filtered_cells = anndata_norm_totalCounts(counts_T, counts_cast_chrX, counts_b6_chrX)
   




######################################################################
# Create table with allelic counts per gene and Xist bin
######################################################################

# Filter for Xist+ cells
adata_XistPos = adata[adata.obs.Xist_pos_neg == 'Xist_pos']
# Convert to df
df_XistPos = adata_XistPos.to_df().join(adata_XistPos.obs)
# Filter for day 3
df_XistPos_day3 = df_XistPos[df_XistPos.index.str.startswith('d3')]
# Filter allelic counts for day3 and Xist+
counts_b6_chrX_filtered_XistPos_day3 = counts_b6_chrX.loc[df_XistPos_day3.index]
counts_cast_chrX_filtered_XistPos_day3 = counts_cast_chrX.loc[df_XistPos_day3.index]
# Melt the dataframes to long format
counts_b6_chrX_filtered_XistPos_day3_long = counts_b6_chrX_filtered_XistPos_day3.reset_index().melt(id_vars='index', var_name='gene', value_name='count_B6')
counts_b6_chrX_filtered_XistPos_day3_long = counts_b6_chrX_filtered_XistPos_day3_long.rename(columns={'index': 'cell'})
counts_cast_chrX_filtered_XistPos_day3_long = counts_cast_chrX_filtered_XistPos_day3.reset_index().melt(id_vars='index', var_name='gene', value_name='count_Cast')
counts_cast_chrX_filtered_XistPos_day3_long = counts_cast_chrX_filtered_XistPos_day3_long.rename(columns={'index': 'cell'})
# Merge the two long format dataframes
allele_counts_merged_day3 = pd.merge(counts_b6_chrX_filtered_XistPos_day3_long, counts_cast_chrX_filtered_XistPos_day3_long, on=['cell', 'gene'])
# Add Xist counts and norm. Xist
Xist_df = pd.DataFrame(adata.obs.total_Xist_counts)
Xist_df['Xist_norm'] = pd.DataFrame(adata.obs.Xist_norm)['Xist_norm']
Xist_df = Xist_df.reset_index().rename(columns={'index': 'cell'})
# Add Xist info to allele df
allele_counts_merged_day3_final = pd.merge(allele_counts_merged_day3, Xist_df, on='cell', how='left')

# Save table
allele_counts_merged_day3_final.to_csv(os.path.join(output_dir, 'AllelicCounts_pacini_day3.txt'), sep='\t', index=False, na_rep='NaN')


