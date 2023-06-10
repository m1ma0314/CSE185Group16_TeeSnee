#If downloading data from 10X genomics...

# We had to add the install path to sys.path to get these imports to work
import os
import sys
sys.path.append(os.environ["HOME"]+"/.local/lib/python3.9/site-packages")

# Import the libraries we installed
import scanpy as sc
import harmonypy
import leidenalg
import anndata as ad

adata = sc.read_10x_mtx('/home/m1ma/teams/16/data/', var_names='gene_symbols', cache=True)
# replace '/home/m1ma/teams/16/data/' with your own directory where you stored your downloaded data

# Inspect the variables (genes/features). you may put it back if you want
# adata.var

# Inspect the observations (cells)
# adata.obs

n_cells_subset = 100  # Specify the desired number of cells in the subset
subset_cells = adata.obs_names[:n_cells_subset]  # Select the first n_cells_subset cell names
adata_subset = adata[subset_cells, :] 

adata_subset.write_csvs("/home/m1ma/teams/16/data/",skip_data=False, sep=',')
# replace '/home/m1ma/teams/16/data/' with your own directory where you want to store your data

# the matrix file is named as X.csv. you could rename it to whatever you want and implement teesnee tool



