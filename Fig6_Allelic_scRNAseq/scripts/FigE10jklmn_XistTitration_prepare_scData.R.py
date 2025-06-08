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
os.chdir(wd)


input_dir = "./input_files/XistTitration_barcodes"
output_dir = "./output_files/"


# Get list of chrX genes
gtf_filename = "./input_files/GENCODE_vM25_plus_Xert_Linx.gtf"
gtf = pd.DataFrame(gtfparse.read_gtf(gtf_filename))
chrX_genes = gtf[gtf[0] == 'chrX'][8].unique() # col 0 is chromosome and col 8 is gene id



######################################################################
# Read-in data
######################################################################

### Read-in count tables
counts_all_rep1 = pd.read_table("./input_files/XistTitration_R1_counts.tsv")
counts_all_rep2 = pd.read_table("./input_files/XistTitration_R2_counts.tsv")
counts_ha1_rep1 = pd.read_table("./input_files/XistTitration_R1_counts_haplo1_chrX.tsv")
counts_ha2_rep1 = pd.read_table("./input_files/XistTitration_R1_counts_haplo2_chrX.tsv")
counts_ha1_rep2 = pd.read_table("./input_files/XistTitration_R2_counts_haplo1_chrX.tsv")
counts_ha2_rep2 = pd.read_table("./input_files/XistTitration_R2_counts_haplo2_chrX.tsv")

### Cell barcodes files
barcode_files_rep1 = {
	'barcodes_titr_Xist_enh_000_diff_4d_R1.txt': 'Xist-RE_dTAG-0',
	'barcodes_titr_Xist_enh_001_diff_4d_R1.txt': 'Xist-RE_dTAG-1',
	'barcodes_titr_Xist_enh_002_diff_4d_R1.txt': 'Xist-RE_dTAG-2',
	'barcodes_titr_Xist_enh_004_diff_4d_R1.txt': 'Xist-RE_dTAG-4',
	'barcodes_titr_Xist_enh_008_diff_4d_R1.txt': 'Xist-RE_dTAG-8',
	'barcodes_titr_Xist_enh_500_diff_4d_R1.txt': 'Xist-RE_dTAG-500',
	'barcodes_titr_Xist_prom_000_diff_4d_R1.txt': 'Xist-PROM_dTAG-0',
	'barcodes_titr_Xist_prom_001_diff_4d_R1.txt': 'Xist-PROM_dTAG-1',
	'barcodes_titr_Xist_prom_002_diff_4d_R1.txt': 'Xist-PROM_dTAG-2',
	'barcodes_titr_Xist_prom_004_diff_4d_R1.txt': 'Xist-PROM_dTAG-4',
	'barcodes_titr_Xist_prom_008_diff_4d_R1.txt': 'Xist-PROM_dTAG-8',
	'barcodes_titr_Xist_prom_500_diff_4d_R1.txt': 'Xist-PROM_dTAG-500'
	}

barcode_files_rep2 = {
	'barcodes_titr_Xist_enh_000_diff_4d_R2.txt': 'Xist-RE_dTAG-0',
	'barcodes_titr_Xist_enh_001_diff_4d_R2.txt': 'Xist-RE_dTAG-1',
	'barcodes_titr_Xist_enh_002_diff_4d_R2.txt': 'Xist-RE_dTAG-2',
	'barcodes_titr_Xist_enh_004_diff_4d_R2.txt': 'Xist-RE_dTAG-4',
	'barcodes_titr_Xist_enh_008_diff_4d_R2.txt': 'Xist-RE_dTAG-8',
	'barcodes_titr_Xist_enh_500_diff_4d_R2.txt': 'Xist-RE_dTAG-500',
	'barcodes_titr_Xist_prom_000_diff_4d_R2.txt': 'Xist-PROM_dTAG-0',
	'barcodes_titr_Xist_prom_001_diff_4d_R2.txt': 'Xist-PROM_dTAG-1',
	'barcodes_titr_Xist_prom_002_diff_4d_R2.txt': 'Xist-PROM_dTAG-2',
	'barcodes_titr_Xist_prom_004_diff_4d_R2.txt': 'Xist-PROM_dTAG-4',
	'barcodes_titr_Xist_prom_008_diff_4d_R2.txt': 'Xist-PROM_dTAG-8',
	'barcodes_titr_Xist_prom_500_diff_4d_R2.txt': 'Xist-PROM_dTAG-500'
	}

# Ensembl IDs of genes of interest to be plotted
genes_of_interest = {
    'Xist': 'ENSMUSG00000086503.4'
    }

# Define what samples will be plotted and in what order
order = ['Xist-RE_dTAG-500','Xist-RE_dTAG-8','Xist-RE_dTAG-4','Xist-RE_dTAG-2','Xist-RE_dTAG-1','Xist-RE_dTAG-0','Xist-PROM_dTAG-500','Xist-PROM_dTAG-8','Xist-PROM_dTAG-4','Xist-PROM_dTAG-2','Xist-PROM_dTAG-1','Xist-PROM_dTAG-0']
sample_prom = ['Xist-PROM_dTAG-500','Xist-PROM_dTAG-8','Xist-PROM_dTAG-4','Xist-PROM_dTAG-2','Xist-PROM_dTAG-1','Xist-PROM_dTAG-0']
sample_re  = ['Xist-RE_dTAG-500','Xist-RE_dTAG-8','Xist-RE_dTAG-4','Xist-RE_dTAG-2','Xist-RE_dTAG-1','Xist-RE_dTAG-0']




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
allele_counts_rep1 = prepareCounts(barcode_files_rep1, 'rep1', counts_ha1_rep1, counts_ha2_rep1)
allele_counts_rep2 = prepareCounts(barcode_files_rep2, 'rep2', counts_ha1_rep2, counts_ha2_rep2)

### Total (non-allelic) counts
# Apply function to filter count tables by cell barcodes and assign the correct samples
total_counts_rep1 = prepareCounts(barcode_files_rep1, 'rep1', counts_all_rep1)
total_counts_rep2 = prepareCounts(barcode_files_rep2, 'rep2', counts_all_rep2)
# Wide format of count table
total_counts_rep1_wide = total_counts_rep1.pivot_table(index=['cell', 'sample'], columns='gene', values='count')
total_counts_rep2_wide = total_counts_rep2.pivot_table(index=['cell', 'sample'], columns='gene', values='count')
# Add row indeces
total_counts_rep1_wide.reset_index(inplace=True)
total_counts_rep2_wide.reset_index(inplace=True)
# Replace NaN with 0
total_counts_rep1_wide = total_counts_rep1_wide.fillna(0)
total_counts_rep2_wide = total_counts_rep2_wide.fillna(0)

print("Total counts prepared")



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
	gene_list = counts_wide.columns.tolist()[2:]  # exclude 'cell' and 'sample'
	sample_list = counts_wide['sample'].tolist()
	cell_list = counts_wide['cell'].tolist()
	X = counts_wide.iloc[:,2:].values

    # Create obs DataFrame to guarantee alignment
	obs_df = pd.DataFrame({
        'sample': sample_list
	}, index=cell_list)

    # Convert counts matrix to anndata object
	adata = ad.AnnData(X, obs=obs_df)
	adata.var_names = gene_list
	adata.obs_names = cell_list

    # Save count data as a layer
	adata.layers["counts"] = adata.X.copy()

    ### Calculate allelic ratios (excluding Xist)
	allele_counts_noXist = allele_counts[~allele_counts['gene'].str.contains(genes_of_interest['Xist'])]
	allele_counts_noXist_perCell = allele_counts_noXist.groupby(['cell', 'sample']).sum().reset_index()
	allele_counts_noXist_perCell = allele_counts_noXist_perCell[['cell', 'sample', 'count_ha1', 'count_ha2']]
	allele_counts_noXist_perCell['allelic_ratio'] = allele_counts_noXist_perCell['count_ha1'] / (
        allele_counts_noXist_perCell['count_ha1'] + allele_counts_noXist_perCell['count_ha2']
	)

    # Merge back into obs using cell as index
	allele_meta = allele_counts_noXist_perCell.set_index('cell')
	adata.obs = adata.obs.join(allele_meta[['allelic_ratio']], how='left')

    ### Extra metadata
	adata.obs['replicate'] = rep
	adata.obs['Xist_raw'] = pd.DataFrame(
	adata[:, adata.var_names == genes_of_interest['Xist']].X.flatten()
	).values
	adata.obs['Xist_pos_neg'] = np.where(adata.obs.Xist_raw > 0, 'Xist_pos', 'Xist_neg')

    ### Mitochondrial QC
	adata.var_names = adata.var_names.astype(str)
	mito_gene_names = sc.queries.mitochondrial_genes("mmusculus", attrname="ensembl_gene_id")
	mito_ids = mito_gene_names["ensembl_gene_id"].tolist()
	adata.var["mt"] = adata.var_names.isin(mito_ids)
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # QC plots
	fig, ax = plt.subplots(figsize=(10, 6))
	sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, show=False, ax=ax)
	fig.savefig(os.path.join(output_dir, "XistTitration_violin_qc_metrics.pdf"), bbox_inches='tight')

    # Filter cells with >3% mito reads
	all_cells = adata.obs_names
	adata = adata[adata.obs.pct_counts_mt < 3, :]
	filtered_cells = np.setdiff1d(all_cells, adata.obs_names)

    ### Normalize counts
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)

    # Add normalized gene expression of genes_of_interest to obs
	for gene in genes_of_interest.keys():
		adata.obs[gene] = pd.DataFrame(
			adata[:, adata.var_names == genes_of_interest[gene]].X.flatten()
		).values

	return adata, filtered_cells

# Apply function to obtain an anndata object with normalized gene expression values
adata_rep1, barcodes_filteredCells_rep1 = anndata_norm_totalCounts(total_counts_rep1_wide, allele_counts_rep1, 'rep1')
print("did first anndata R1")
adata_rep2, barcodes_filteredCells_rep2 = anndata_norm_totalCounts(total_counts_rep2_wide, allele_counts_rep2, 'rep2')
print("did first anndata R2")




######################################################################
# Filter out XO cells and combine both replicates; filter for Xist+ cells
######################################################################

### Filter out XO cells from anndata object
adata_rep1_noXO = adata_rep1[~((adata_rep1.obs.Xist < 0.5) & ((adata_rep1.obs.allelic_ratio < 0.15) | (adata_rep1.obs.allelic_ratio > 0.85)))]
adata_rep2_noXO = adata_rep2[~((adata_rep2.obs.Xist < 0.5) & ((adata_rep2.obs.allelic_ratio < 0.15) | (adata_rep2.obs.allelic_ratio > 0.85)))]

##Print the number of XO cells by sample
adata_rep1_XO = adata_rep1[((adata_rep1.obs.Xist < 0.5) & ((adata_rep1.obs.allelic_ratio < 0.15) | (adata_rep1.obs.allelic_ratio > 0.85)))]
adata_rep2_XO = adata_rep2[((adata_rep2.obs.Xist < 0.5) & ((adata_rep2.obs.allelic_ratio < 0.15) | (adata_rep2.obs.allelic_ratio > 0.85)))]
print(adata_rep1_XO.obs['sample'].value_counts().sort_index())
print(adata_rep2_XO.obs['sample'].value_counts().sort_index())


# Get the barcodes of the XO cells which get filtered out
barcodes_filteredCells_XO_rep1 = np.array(adata_rep1[(adata_rep1.obs.Xist < 0.5) & ((adata_rep1.obs.allelic_ratio < 0.15) | (adata_rep1.obs.allelic_ratio > 0.85))].obs_names)
barcodes_filteredCells_XO_rep2 = np.array(adata_rep2[(adata_rep2.obs.Xist < 0.5) & ((adata_rep2.obs.allelic_ratio < 0.15) | (adata_rep2.obs.allelic_ratio > 0.85))].obs_names)

print("XO filtered")

### Combine the two replicates
adata_combined = ad.concat([adata_rep1_noXO, adata_rep2_noXO], axis=0, join='outer')
# Replace NaN values with 0 in the .X matrix
adata_combined.X = np.nan_to_num(adata_combined.X)

print("adata combined")

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

print("Xist+ filtered")



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
plot_pca(adata_combined, 'sample', None, None, filename= os.path.join(output_dir, "FigE10j_XistTitration_sample_PCA.pdf"))
plot_pca(adata_combined, 'replicate', None, None, filename= os.path.join(output_dir, "FigE10k_XistTitration_rep_PCA.pdf"))

def plot_pca_faceted(adata, palette=None, colormap=None, filename=None):
	
	# HVGs and PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=500, batch_key="sample")
    sc.tl.pca(adata)

    # Create figure
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # PROM samples
    adata_prom = adata[adata.obs['sample'].isin(sample_prom)].copy()
    sc.pl.pca(adata_prom, color='sample', ax=axs[0], show=False, palette=palette, color_map=colormap, size=40)
    axs[0].set_title('Xist-PROM Samples')

    # RE samples
    adata_re = adata[adata.obs['sample'].isin(sample_re)].copy()
    sc.pl.pca(adata_re, color='sample', ax=axs[1], show=False, palette=palette, color_map=colormap, size=40)
    axs[1].set_title('Xist-RE Samples')

    plt.tight_layout()
    if filename:
        fig.savefig(os.path.join(output_dir, filename), bbox_inches='tight')
    plt.show()


plot_pca_faceted(adata_combined, palette="Set2", filename= os.path.join(output_dir, "FigE10l_XistTitration_dTAG_PCA.pdf"))



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
    # Allelic QC plots
	features = ["n_genes_by_counts", "total_counts"]
	fig, axs = plt.subplots(1, len(features), figsize=(6 * len(features), 5))

	for i, feature in enumerate(features):
		sc.pl.violin(
			adata,
			keys=feature,
			jitter=0.4,
			ax=axs[i],
			show=False
			)
	fig.tight_layout()
	fig.savefig(os.path.join(output_dir, "FigE10mn_XistTitration_allelic_qc_metrics.pdf"), bbox_inches='tight')
    
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
df.to_csv(os.path.join(output_dir, 'XistTitration_CountTable_plusMetadata_totalCounts.txt'), sep='\t', index=True, na_rep='NaN')

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
allele_counts_XistPos_concat_plusXist.to_csv(os.path.join(output_dir, 'XistTitration_AllelicCounts_plusTotalXist.txt'), sep='\t', index=True, na_rep='NaN')
