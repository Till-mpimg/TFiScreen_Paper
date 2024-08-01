#!/usr/bin/env python3

import pickle
import pandas as pd
import sklearn
import numpy as np
import scanpy as sc
import scipy.sparse
import anndata
import matplotlib
import matplotlib.pylab as plt
import os 
import sys
import warnings
import seaborn as sns
import random as rn 
from sklearn.decomposition import PCA
import copy
import scATAcat

# set seed for reproducibility
sd = 1234
np.random.seed(sd)
rn.seed(sd)
%env PYTHONHASHSEED=0

#Set work directory
wd = sys.argv[1]
os.chdir(wd)

#Detail directories
results_dir = "./output_files/"
data_dir = "./data/"
input_dir = "./input_files/"

scRNA_count = pd.read_table(data_dir + "mohammed_2017_counts.txt", sep="\t")

#exclude genes expressed in less than three cells
scRNA_count_filtered = scRNA_count[scRNA_count.iloc[:,3:].ne(0).sum(1).gt(2)]

# filter autosomal counts
scRNA_count_filtered_autosomal = scRNA_count_filtered[~scRNA_count_filtered.seqnames.isin(['chrM','chrX','chrY'])]
scRNA_count_filtered_autosomal_matrix = scRNA_count_filtered_autosomal.drop(['seqnames', 'gene_id'], axis=1)
scRNA_count_filtered_autosomal_matrix.drop(['gene_name'], axis=1, inplace=True)
scRNA_count_filtered_autosomal_matrix[~scRNA_count_filtered_autosomal_matrix.index.duplicated(keep='first')]


# initialize AnnData matrix
genes_df =pd.DataFrame({"genes":scRNA_count_filtered_autosomal_matrix.index } , index=scRNA_count_filtered_autosomal_matrix.index)
cells_df =pd.DataFrame({"cells":scRNA_count_filtered_autosomal_matrix.columns },index=scRNA_count_filtered_autosomal_matrix.columns)
scRNA_count_uniq_csr= scipy.sparse.csr_matrix(scRNA_count_filtered_autosomal_matrix.T.astype(pd.SparseDtype("int64",0)).sparse.to_coo())
sc_adata= anndata.AnnData(scRNA_count_uniq_csr, var=genes_df, obs=cells_df)

# preprocessing is not necessary because the input data is already QC fitered
# normalize and scale the data
sc.pp.normalize_total(sc_adata, target_sum=1e6)
sc_adata.X= np.array(np.log2(sc_adata.X.todense()+1))

#PCA embedding 
#Get timepoints from cell name
sc.tl.pca(sc_adata, svd_solver='arpack')
sc_adata.obs["timepoint_sc"] = [i.split("_")[0] for i in sc_adata.obs["cells"]]

#Plots PCA colored by timepoint
sc.pl.pca(sc_adata, color='timepoint_sc')

#For each timepoint we want to exclude extraembryonic cells
#Take all cells from E3.5 (before separation)


# E4.5 epi markers Nanog, Esrrb, Fgf4 (epiblast)
## At E4.5, a clear separation of cells into the epiblast and PrE is 
# observed and characterized by exclusive expression of known markers,
# such as Nanog, Esrrb, Fgf4 (epiblast) 
# and Gata6, Pdgfra, and Sox17 (PrE) 
sc.pl.pca(sc_adata ,color=["timepoint_sc",'Nanog','Esrrb', 'Fgf4'])

sc_adata_E4dot5= sc_adata[sc_adata.obs.timepoint_sc.isin(['E4.5']) ]
E4dot5_epi_cells = sc_adata_E4dot5[(sc_adata_E4dot5[: , 'Nanog'].X > 3) | (sc_adata_E4dot5[: , 'Esrrb'].X > 3) | (sc_adata_E4dot5[: , 'Fgf4'].X > 3), :].obs 
sc.pl.pca(sc_adata_E4dot5, color=["timepoint_sc",'Nanog','Esrrb', 'Fgf4','Zfp57', "Pou5f1", "Sox17", "Gata4", "Gata6", "T"])

#Thresholds are retrieved from visual inspection of gene expression histograms
sns.displot(sc_adata_E4dot5[:, 'Fgf4'].X,kde=True, bins=15)
E4dot5_epi_cells = sc_adata_E4dot5[(sc_adata_E4dot5[: , 'Nanog'].X > 3) | (sc_adata_E4dot5[: , 'Esrrb'].X > 3) | (sc_adata_E4dot5[: , 'Fgf4'].X > 3), :].obs 
sc_adata.obs['E4dot5_epi'] =  np.where(sc_adata.obs.index.isin(E4dot5_epi_cells.index), 'E4.5 epiblasts', 'other cells')
sc_adata.obs['E4dot5_epi']= sc_adata.obs['E4dot5_epi'].astype('category')


# The E5.5 epiblast cells cluster separately from E4.5 epiblast cells and 
#possess reduced Nanog expression, while gaining primed pluripotency markers 
# such as Pou3f1. Cells within the E5.5 epiblast do not show any apparent substructure
sc_adata_E5dot5= sc_adata[sc_adata.obs.timepoint_sc.isin(['E5.5']) ]

# E5.5 epi markers
sc.pl.pca(sc_adata_E5dot5, color=["timepoint_sc",'Gm5529','Psme2', "Dnmt3a", "Hmga1","Pou3f1",'Nanog','Esrrb', 'Fgf4', "Pou5f1", "Sox17", "Gata4", "Gata6", "Gata2", "Gata3"])

# Filter epi cells
E5dot5_epi_cells = sc_adata_E5dot5[(sc_adata_E5dot5[: , 'Dnmt3a'].X > 7) & (sc_adata_E5dot5[: , 'Pou3f1'].X > 3) , :].obs 
sc_adata_E5dot5.obs['E5dot5_epi'] =  np.where(sc_adata_E5dot5.obs.index.isin(E5dot5_epi_cells.index), 'E5.5 epiblasts', 'other E5.5 cells')
sc_adata_E5dot5.obs['E5dot5_epi']= sc_adata_E5dot5.obs['E5dot5_epi'].astype('category')
sc_adata.obs['E5dot5_epi'] =  np.where(sc_adata.obs.index.isin(E5dot5_epi_cells.index), 'E5.5 epiblasts', 'other cells')
sc_adata.obs['E5dot5_epi']= sc_adata.obs['E5dot5_epi'].astype('category')


#E6.5/E6.75
#Remove visceral endoderm, include primitive streak
sc_adata_E6dot5= sc_adata[sc_adata.obs.timepoint_sc.isin(['E6.5', 'E6.75']) ]

# E6.5 epi1 markers : Gm6984 Zic5
sc.pl.pca(sc_adata_E6dot5, color=["timepoint_sc",'Zic5',"Pou5f1", "Sox17", "Gata4", "Gata6", "Ttr", "Ctsh", "T"])



# E6.5 Viseral endoderm
sc.pl.pca(sc_adata_E6dot5, color=["timepoint_sc","Lhx1","Peg3","Eomes","Creb3l3","Maf"])

# E6.5 primitive streak
sc.pl.pca(sc_adata_E6dot5, color=["timepoint_sc","Snai1", "Lef1","Evx1", "Mesp1", "T"])

#Check histograms to set threshold
sns.displot(sc_adata_E6dot5[:, 'Lhx1'].X,kde=True, bins=15)
E6dot5_epi_cells = sc_adata_E6dot5[~((sc_adata_E6dot5[: , 'Maf'].X > 2) |(sc_adata_E6dot5[: , 'Creb3l3'].X > 2) |(sc_adata_E6dot5[: , 'Lhx1'].X > 8))  ,].obs


sc_adata_E6dot5.obs['E6dot5_epi'] =  np.where(sc_adata_E6dot5.obs.index.isin(E6dot5_epi_cells.index), 'E6.5 epiblasts', 'other E6.5 cells')
sc_adata_E6dot5.obs['E6dot5_epi']= sc_adata_E6dot5.obs['E6dot5_epi'].astype('category')

sc_adata.obs['E6dot5_epi'] =  np.where(sc_adata.obs.index.isin(E6dot5_epi_cells.index), 'E6.5 epiblasts', 'other E6.5 cells')
sc_adata.obs['E6dot5_epi']= sc_adata.obs['E6dot5_epi'].astype('category')


#Combine epi cells into new sc object
epi_cells = sc_adata.obs[sc_adata.obs['timepoint_sc']==("E3.5")].cells +E4dot5_epi_cells.cells + E5dot5_epi_cells.cells +E6dot5_epi_cells.cells
epi_cells = list(epi_cells.index)
sc_adata_epi = sc_adata[epi_cells].copy()

# remane timepoint E6.75 s 6.5
rename_timepoint_dict = {"E6.75":"E6.5"}
sc_adata.obs = sc_adata.obs.replace({"timepoint_sc": rename_timepoint_dict})

# remane timepoint E6.75 s 6.5
rename_timepoint_dict = {"E6.75":"E6.5"}
sc_adata_epi.obs = sc_adata_epi.obs.replace({"timepoint_sc": rename_timepoint_dict})

#Integrate with bulk data
bulk_rnaseq_counts = pd.read_table(input_dir + "Counts_RNA_timecourse.txt", sep="\t")
gene_list = list(bulk_rnaseq_counts["gene"])


#make gene names unique
unique_gene_names =  [v + "_" + str(gene_list[:i].count(v) + 1) if gene_list.count(v) > 1 else v for i, v in enumerate(gene_list)]
duplicating_genes = set([x for x in gene_list if gene_list[0:1000].count(x) > 1])

# keep the gene sym ensID as a df:
gene_id_mapping_df = pd.DataFrame({'org_gene_sym' : bulk_rnaseq_counts["gene"], "ENS_ID" : bulk_rnaseq_counts["gene_id"], "unique_gene_sym": unique_gene_names})
bulk_rnaseq_counts.index = unique_gene_names
bulk_rnaseq_counts.drop(["gene", "gene_id"], axis=1,inplace=True)     

#Use scATAcat to built AnnData
bulk_adata = scATAcat.generate_bulk_sparse_AnnData(bulk_rnaseq_counts)
bulk_adata.var.columns = ["genes"]

#Overlap bulk and sc features
#First remove features without count and filter for common features
sc_bulk_common_vars = scATAcat.overlap_vars(sc_adata_epi, bulk_adata)
sc_epi_commonFeatures_adata = sc_adata_epi[:, sc_adata_epi.var.genes.isin(sc_bulk_common_vars)]
bulk_commonFeatures_adata = scATAcat.subset_adata_vars(bulk_adata, vars_list=sc_bulk_common_vars, copy_=True)

#Subset to differential genes
diffGenes = pd.read_csv(data_dir + "alldiffGenes_timeseries.csv", index_col=0)

common_differential_vars = list(set(list(sc_bulk_common_vars)) & set(list(diffGenes['x'])))
bulk_commonDiffFeatures_adata = scATAcat.subset_adata_vars(bulk_commonFeatures_adata, vars_list=common_differential_vars, copy_=True)


sc_epi_commonDiffFeatures_adata = scATAcat.subset_adata_vars(sc_epi_commonFeatures_adata, vars_list=common_differential_vars, copy_=True)

# add cell_type  as obs
sc_epi_commonDiffFeatures_adata.obs["cell_type"] = [r.split('_')[0] for r in sc_epi_commonDiffFeatures_adata.obs.cells]


#Prepare anndata
sc_epi_commonFeatures_adata_raw = anndata.AnnData(scRNA_count_uniq_csr, var=genes_df, obs=cells_df)
sc_epi_commonFeatures_adata_raw = sc_epi_commonFeatures_adata_raw[sc_epi_commonFeatures_adata.obs.cells].copy()
sc_epi_commonFeatures_adata_raw = sc_epi_commonFeatures_adata_raw[:,sc_epi_commonFeatures_adata_raw.var.genes.isin(sc_epi_commonFeatures_adata.var.genes)].copy()


# add cell_type  as obs
sc_epi_commonFeatures_adata_raw.obs["timepoint_sc"] = [r.split('_')[0] for r in sc_epi_commonFeatures_adata_raw.obs.cells]
rename_timepoint_dict = {"E6.75":"E6.5"}
sc_epi_commonFeatures_adata_raw.obs = sc_epi_commonFeatures_adata_raw.obs.replace({"timepoint_sc": rename_timepoint_dict})

#Pseudobulk by timepoint
pseudobulk_epi_commonFeatures_adata = scATAcat.generate_bulk_sparse_AnnData(scATAcat.get_pseudobulk_matrix_ext(adata_to_subset=sc_epi_commonFeatures_adata_raw, adata_to_get_clusters=sc_epi_commonFeatures_adata_raw, cluster_key="timepoint_sc"  ,method = 'sum'))
scATAcat.preprocessing_libsize_norm_log2(pseudobulk_epi_commonFeatures_adata)
scATAcat.preprocessing_libsize_norm_log2(bulk_commonFeatures_adata)

# Subset to diff regions and z normalize
bulk_commonDiffFeatures_adata = scATAcat.subset_adata_vars(bulk_commonFeatures_adata, vars_list=common_differential_vars, copy_=True)
pseudobulk_epi_commonDiffFeatures_adata = scATAcat.subset_adata_vars(pseudobulk_epi_commonFeatures_adata, vars_list=common_differential_vars, copy_=True)
scATAcat.preprocessing_standardization(bulk_commonDiffFeatures_adata, input_layer_key="libsize_norm_log2", zero_center=True)
scATAcat.preprocessing_standardization(pseudobulk_epi_commonDiffFeatures_adata, input_layer_key="libsize_norm_log2", zero_center=False, output_layer_key= "libsize_norm_log2_bulk_scaled_diff",
                              std_key= None,  mean_key=None, std_ = bulk_commonDiffFeatures_adata.var["feature_std"], mean_= bulk_commonDiffFeatures_adata.var["feature_mean"])

#Run projection of both files
def projection_match_colors_sexed_2D(bulk_adata, pseudobulk_adata, bulk_layer_key = "libsize_norm_log2_std", pseudobulk_layer_key="libsize_norm_log2_bulk_scaled", color_key="clustering_color", label_font_size =18):
    num_of_bulk_samples = bulk_adata.n_obs
    if num_of_bulk_samples >30:
        n_comp = 30
    else:
        n_comp = num_of_bulk_samples
    pca_bulk = PCA(n_components=n_comp)
    pca_bulk_train = pca_bulk.fit_transform(bulk_adata.layers[bulk_layer_key])
    PCs_pseudobulk_projection = pca_bulk.transform(pseudobulk_adata.layers[pseudobulk_layer_key])


    trained_bulk_pca_df = pd.DataFrame(data = pca_bulk_train
                 , columns = ["principal component " + str(i) for i in range(n_comp)])
    trained_bulk_pca_df["targets"] = bulk_commonDiffFeatures_adata.obs.index
    trained_bulk_pca_df

    projected_pseudobulk_pca_df = pd.DataFrame(data = PCs_pseudobulk_projection
                 , columns = ["principal component " + str(i) for i in range(n_comp)])
    projected_pseudobulk_pca_df["targets"] = pseudobulk_adata.obs.index
    projected_pseudobulk_pca_df

    trained_bulk_pca_df_w_labels = trained_bulk_pca_df.copy()
    trained_bulk_pca_df_w_labels["cell_type"] = list((['_'.join([i.split("_")[2], i.split("_")[0]]) for i in trained_bulk_pca_df_w_labels['targets']]))
    #trained_bulk_pca_df_w_labels['targets'].apply(lambda r: '_'.join(r.split('_')[:-1]))
    my_color=pd.Series(pd.Categorical(trained_bulk_pca_df_w_labels["cell_type"])).cat.codes
    trained_bulk_pca_df_w_labels["color_id"] = my_color
    ## plot
    fig = plt.figure()
    ax = fig.add_subplot(111)

    PC1 = trained_bulk_pca_df_w_labels['principal component 0'].values
    PC2 = trained_bulk_pca_df_w_labels['principal component 1'].values

    CELLTYPES = trained_bulk_pca_df_w_labels["cell_type"].values
    CELLTYPES_ = np.unique(CELLTYPES)
    num_bulk_cell_types= len(CELLTYPES_)
    # colors for bulk samples 

    COLORS = ["#f9de21",
    "#f7e225",
    "#fccd25",
    "#feb72d",
    "#fca338",
    "#f79044",
    "#f07f4f",
    "#e76e5b",
    "#dd5e66",
    "#d14e72",
    "#c5407e",
    "#b6308b",
    "#a72197",
    "#9511a1",
    "#8305a7",
    "#6e00a8",
    "#5901a5",
    "#43039e",
    "#2c0594",
    "#0d0887"]

    for cell_type1, color in zip(CELLTYPES_, COLORS):
        idxs = np.where(CELLTYPES == cell_type1)
        # No legend will be generated if we don't pass label=species
        ax.scatter(
            PC1[idxs,], PC2[idxs,], label=cell_type1,
            s=200, color=color, alpha=0.7, marker = "v"
        )
    m= np.array([list(PC1),list(PC2)])
    for i in range(len(m[0])): #plot each point + it's index as text above
        ax.text(m[0,i],m[1,i],  '%s' % ('  ' +CELLTYPES[i]), fontsize=label_font_size)
    
    BULK_CELLTYPES = projected_pseudobulk_pca_df["targets"].values
    CELLTYPES_ = np.unique(BULK_CELLTYPES)
    PC1_ = projected_pseudobulk_pca_df['principal component 0'].values
    PC2_ = projected_pseudobulk_pca_df['principal component 1'].values
    # # get the color code frompseudobulk matrix:
    COLORS_ =pseudobulk_adata.uns[color_key]
    #["r", "b","g", "tab:pink", "tab:cyan" ] #'#A05529'
    for bulk_cell_type1, color in zip(CELLTYPES_, COLORS_):
        #print(CELLTYPES_)
        #print(COLORS_)
        idxs = np.where(BULK_CELLTYPES == bulk_cell_type1)
        # No legend will be generated if we don't pass label=species
        ax.scatter(
            PC1_[idxs,], PC2_[idxs,], label=bulk_cell_type1,
            s=150, color=color, alpha=0.9)
        #ax.scatter(trained_pseudobulk_pca_df['principal component 1'], trained_pseudobulk_pca_df['principal component 2'], trained_pseudobulk_pca_df['principal component 3'], c=[4,5,6,7,8] ,cmap="Set2_r", s=60, marker = 'v')

    # put labels on the plot
    m= np.array([list(PC1_),list(PC2_)])
    for i in range(len(m[0])): #plot each point + it's index as text above
        ax.text(m[0,i],m[1,i],  '%s' % ('  ' +BULK_CELLTYPES[i]), fontsize=label_font_size)


    #fig
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax.legend( ncol=2,handleheight=2.4, labelspacing=0.05, bbox_to_anchor=(0.9, 0.5, 0.8, 0.5))
    figure = plt.gcf() # get current figure
    figure.set_size_inches(20 ,30)
    return fig, trained_bulk_pca_df_w_labels,projected_pseudobulk_pca_df

res_matching_sexed = projection_match_colors_sexed(bulk_commonDiffFeatures_adata, pseudobulk_epi_commonDiffFeatures_adata, bulk_layer_key = "libsize_norm_log2_std", pseudobulk_layer_key="libsize_norm_log2_bulk_scaled_diff" ,color_key= 'leiden_colors')

#Output pca info
pca_bulk = res_matching_sexed[1]
pca_pseudobulk = res_matching_sexed[2]
pca_bulk.to_csv(results_dir + "bulk_pca.tsv", sep="\t")
pca_pseudobulk.to_csv(results_dir + "pseudobulk_pca.tsv", sep="\t")


#Output distance matricx
centroid_heatmap =  scATAcat.plot_pca_dist_cent_heatmap(res_matching_sexed[1],res_matching_sexed[2])
distance_matrix = scATAcat.get_pseudobulk_to_prototype_distance(centroid_heatmap[1], pbulk_to_prototype=False).T
distance_matrix.to_csv(results_dir + "distance_matrix.tsv", sep="\t")

#Write anndata outputs
sc_epi_commonDiffFeatures_adata.write_h5ad(filename= results_dir + "sc_epi_commonDiffFeatures_adata.h5ad")
sc_epi_commonFeatures_adata.write_h5ad(filename=results_dir + "sc_epi_commonFeatures_adata.h5ad")
sc_adata.write_h5ad(filename=results_dir + "sc_adata.h5ad")
bulk_commonDiffFeatures_adata.write_h5ad(filename=results_dir + "bulk_commonDiffFeatures_adata.h5ad")
pseudobulk_epi_commonDiffFeatures_adata.write_h5ad(filename=results_dir + "pseudobulk_epi_commonDiffFeatures_adata.h5ad")










