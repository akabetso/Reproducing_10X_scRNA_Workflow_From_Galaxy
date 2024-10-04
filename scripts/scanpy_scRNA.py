# Import libraries
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

# Define directories
matrix_dir = "../results/DropletUtils"
save_dir = "../results/scanpy"

# Create directory
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Configure Scanpy library for logging verbosity
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.figdir = save_dir

# tDefine the h5ad file to contain the analysis results
results_file = "Input_3K_PBMC.h5ad"  

# Read the matrix directory and parse it to adata object
adata = sc.read_10x_mtx(
    matrix_dir,  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)

# Filter genes found in less than 3 cells
sc.pp.filter_genes(adata, min_cells=3)
 
# Annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# Create a violin plot by genes counts and mt counts using jitter to add random noise for better visibility.
# Violin plots provides insights into quality control by observing the distrubutions accross accross all cells.
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save = "violin_qc_plot.png",
    show = False
)

# Observe the mt counts to the total counts and genes to total counts
# Helps detect outliers and guide downstream analysis decisions.
# Totoal counts is the number of unique RNA
# mt is the transcripts comming from the mitochrondrial genes. Higher mt expresses poor quallity cells or means the cells are stressed or dying
# genes by total counts is the number of unique genes expressed in a cell. lower numbers suggest poor quality cells
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax = ax[0], show = False)
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", ax = ax[1], show = False)
save_path = os.path.join(save_dir, "scatterplot_mt_genes_total_counts.png")
plt.savefig(save_path)
plt.close()

# From the QC observation, we filter cells below 200 and above 2500. We equally cut of mt above 5%.
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
adata = adata[adata.obs.n_genes_by_counts > 200, :]

# Normalization to account for differences in sequencing depth between cells. ensures the total number of counts for each cell is comparable.
# The total counts for each cell should equal 10000. provides balance between maintaining signal and reducing noise.
sc.pp.normalize_total(adata, target_sum=1e4, key_added='norm')

# Perform a logarithmic transformation
# reduce skewness making data more easy to visualise, stabilize variance.
sc.pp.log1p(adata)

#extract and plot highly variable genes with minimum and maximum threshold for the mean dispersion of genes, with a minimum dispersion.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor='seurat')
fig, ax = plt.subplots(1, 2, figsize = (12, 6))
sc.pl.highly_variable_genes(adata, show = False)
save_path = os.path.join(save_dir, "highly_variable_genes.png")
plt.savefig(save_path)
plt.close()

# freeze the current state
adata.raw = adata

# Perform a linear regression to remove bias effects from total counts per cell and mt counts
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

# Scale the data object to unit variance and zero mean. clips value exceeding a standard dev of 10
# The variance accross cells is 1 and the mean expression is 0 which gives equal weith in downstream analysis
# Ensures highly expressed genes do not dominate
sc.pp.scale(adata, max_value=10)

# Perform dimension reduction with PCA
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
# Visualizing PCA
fig, ax = plt.subplots(1, 2, figsize = (12, 6))
sc.pl.pca(adata, components=['1,2'], color="CST3", ax = ax[0], show = False)
sc.pl.pca(adata, components=['2,3'], color="CST3", ncols=2, ax = ax[1], show = False)
save_path = os.path.join(save_dir, "PCA_plot.png")
plt.savefig(save_path)
plt.close()

# Compute the number of PCs to keep using an elbow plot.
# We used 50 PCs for the PCA analysis representing a robust compression of the dataset
sc.pl.pca_variance_ratio(adata, log=True, save = "elbow_plot.png", show = False)
# Compute a neighborhood graph of observations
# Use 10 manifold approximation or neighbors for a KNN graph and 10 PCs as defined by the elbow plot
# Important steps before clustering, gives higher weith to similar expresion profiles.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, method='umap', metric='euclidean')

# Embed and plot the neighborhood placing similar cells together in low dimensional space
sc.tl.umap(adata)

# Clustering the cells within a neighborhood graph
# 0.4 resolution is recomended by Galaxie
sc.tl.louvain(adata, resolution=0.45, flavor='vtraag')

# Save current state of clustered adata object
adata.write("_umap_clustering.h5ad")

# Visualize the clusters
sc.pl.umap(adata, color=["louvain", "CST3", "NKG7", "PPBP"], save = "_louvain_clustering.png", show = False)

# We use t-test to indefity genes that drive seperation between clusters and rank them
sc.tl.rank_genes_groups(adata, groupby='louvain', use_raw=True, method='t-test', n_genes=100, corr_method='benjamini-hochberg')
sc.settings.verbosity = 2  # reduce the verbosity
sc.pl.rank_genes_groups(adata, n_genes=20, ncols=3, sharey=True, save = "_t_test.png", show=False)
# List and compare top 5 elements of the ranked genes
marker_genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)
fig, ax = plt.subplots(figsize=(6, 1))
ax.axis("tight")
ax.axis("off")
table = ax.table(cellText = marker_genes.values, colLabels = marker_genes.columns, cellLoc = "center", loc = "center")
save_path = os.path.join(save_dir, "t_test_ranked_genes_accross_cells.png")
plt.savefig(save_path)
plt.close

# Save current state of adata object
adata.write("_umap_clustering_t_test.h5ad")

# Read the previously clustered data and perform Wilcoxon rank sum test
# This is another widely used method for pairwise comparisons between groups of observations. 
# Directly assesses seperation between the expression distributions of different clusters.
# Will compare results to previous t-test
adata = sc.read_h5ad("_umap_clustering.h5ad")
sc.tl.rank_genes_groups(adata, groupby='louvain', use_raw=True, method='wilcoxon', n_genes=100, corr_method='benjamini-hochberg')
sc.pl.rank_genes_groups(adata, n_genes=20, ncols=3, sharey=True, save = "_wilcoxon.png", show=False)
# List and compare top 5 elements of the wilcoon ranked genes
marker_genes = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5)
fig, ax = plt.subplots(figsize=(5, 1))
ax.axis("tight")
ax.axis("off")
table = ax.table(cellText = marker_genes.values, colLabels = marker_genes.columns, cellLoc = "center", loc = "center")
save_path = os.path.join(save_dir, "wilcoxon_ranked_genes_accross_cells.png")
plt.savefig(save_path)
plt.close

adata.write("_umap_clustering_wilcoxon.h5ad") # Write current state to h5ad file

# Compare different expresion level of define genes accross cells
sc.pl.violin(adata, keys=['CST3', 'NKG7', 'PPBP'], groupby='louvain', use_raw=True, save = "_expresion_level.png", show = False)
marker_genes = ['LDHB', 'CD74', 'LYZ', 'CCL5', 'LST1', 'NKG7', 'HLA-DPA1', 'PF4']
sc.pl.stacked_violin(adata, marker_genes, groupby="louvain", figsize = (8, 8), swap_axes = True, stripplot = False, save = "_mkg_expression_level.png", show = False)
sc.pl.dotplot(adata, marker_genes, groupby="louvain", figsize = ( 10, 6), swap_axes = True, save = "_mkg_expresion_lev.png", show = False)
# Display the mean expression of the marker genes for the cells on the neighborhood graph
sc.pl.umap(adata, 
           color=['louvain', 'LDHB', 'CD74', 'LYZ', 'CCL5', 'LST1', 'NKG7', 'HLA-DPA1', 'PF4'], 
           use_raw=True, 
           ncols=3, save = "_mkg_expression_lev.png", show = False)  # Number of panels per row
# We can equally display the top 20 marker genes in different cells using heatmap
# Plot ranking of genes as a heatmap
sc.pl.rank_genes_groups_heatmap(adata, 
                                 n_genes=20, 
                                 use_raw=True, 
                                 dendrogram=True, save = "_mkg_expression_lev.png", show = False)

# Rank genes using the Wilcoxon Rank-Sum test for groups 0 vs 1
sc.tl.rank_genes_groups(adata, 
                         groupby='louvain', 
                         use_raw=True, 
                         groups=['0'],  
                         reference='1',  
                         n_genes=100, 
                         method='wilcoxon', p_val='Benjamini_Hochberg')  
sc.pl.rank_genes_groups(adata, 
                         n_genes=20, 
                         sharey=False, show = False)  
# Generate the violin plot for the top 10 ranked genes
plt.figure(figsize=(10, 6)) 
sc.pl.rank_genes_groups_violin(
    adata, 
    n_genes=10,  
    use_raw=True, show = False
)
save_path = os.path.join(save_dir, "violin_ranked_genes_0_vs_1.png")
plt.savefig(save_path)
plt.close()

# Perform cell type annotation or annotating each cluster
# Define the mapping of old categories to new categories
old_categories = adata.obs['louvain'].unique()
print(old_categories)
#new_categories = [0, 1, 2, 3, 4, 5, 6, 7]
new_categories = ['CD4+ T', 'B', 'CD14+', 'NK', 'CD8+ T', 'FCGR3A+',  'Dendritic', 'Megakaryocytes']
# Create a mapping dictionary
category_mapping = {old: new for old, new in zip(old_categories, new_categories)}
print(category_mapping)
# Rename categories in the observation annotation
adata.obs['louvain'] = adata.obs['louvain'].map(category_mapping)
# Visualize cluster annotation with umap
sc.pl.umap(
    adata, 
    color='louvain',
    use_raw=True,   
    frameon=False,    
    legend_loc='on data', 
    save = "_annotaion.png",
    show = False 
)

# Thank you for reading through.
# This pipeline was presented in Galaxy
# The workflow used here was produced by carefuly adapting from the documentations of the tools mentioned in Galaxie to reproduce their results
# Presents a horlistic understanding of RNA seq pipelines.



