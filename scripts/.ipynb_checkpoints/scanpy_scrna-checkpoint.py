import numpy as np
import pandas as pd
from scipy.io import mmread
import anndata as ad
import scanpy as sc


# Directories
matrix_dir = "../results/DropletUtils/matrix.mtx"  # Path to your .mtx file
features_dir = "../results/DropletUtils/features.tsv"
barcodes_dir = "../results/DropletUtils/barcodes.tsv"

# Load the sparse matrix from the .mtx file
matrix = mmread(matrix_dir).tocsc()  # Convert to Compressed Sparse Column format

# Load gene names
genes = pd.read_csv(features_dir, sep='\t', header=None, names=['gene_id', 'gene_name'])
gene_names = genes['gene_name'].values

# Load barcodes
barcodes = pd.read_csv(barcodes_dir, sep='\t', header=None, names=['barcode'])
barcode_names = barcodes['barcode'].values

# Create the AnnData object
adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=barcode_names), var=pd.DataFrame(index=gene_names))

# Save the AnnData object to a file
adata.write('3k_PBMC_Input.h5ad')

print("AnnData object created and saved successfully.")