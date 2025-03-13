## MuTrans untreated: GRNB5 data prep
## February 2025

import sys 
import os 
import pandas as pd 
import numpy as np 
import anndata as ad 
import scanpy as sc 
import matplotlib.pyplot as plt 
import hdf5plugin

current_dir = os.getcwd() 
script_dir = "/data/rds/DMP/DUDMP/PAEDCANC/echen/MuTrans-release-main/Example" 
datadir = "data_2025_ut/" 
plotdir = "../plots/2025_ut/" 
filename = "grnb5_seurat.h5ad" 
file_path = os.path.join(datadir, filename) 
sample_name = filename.replace("_seurat.h5ad", "") 
os.chdir(script_dir) 
sc.settings.set_figure_params(dpi=100, frameon=False) 
print(f"Current directory: {current_dir}", flush=True) 
print(f"Changed directory: {os.getcwd()}", flush=True) 
print(f"File directory: {file_path}", flush=True) 

print(f"Processing {filename}...", flush=True)         
adata = sc.read_h5ad(file_path) 
adata.obs = adata.obs[['Sample', 'Barcode', 'batch', 'description', 
                       'Condition', 'Rec', 'seurat_clusters.0.4', 
                       'MES.Sig', 'ADRN.Sig', 'AMT.Sig', 
                       'AMT.score', 'AMT.state']] 
adata.obsm = []         
df = pd.DataFrame.sparse.from_spmatrix(adata.X) 
print(f"Raw counts added for {sample_name}", flush=True) 
print(df.head(10), flush=True)         

print(f"Beginning preprocessing of {sample_name} AnnData object", flush=True) 
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars = ['mt'], percent_top = None, log1p = False, inplace = True)         
print(f"Normalising {sample_name} data", flush=True) 
sc.pp.normalize_total(adata, target_sum = 1e4) 
sc.pp.log1p(adata) 
sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) 
sc.pl.highly_variable_genes(adata, show = False) 
print(f"Save fig to {plotdir}{sample_name}_hvgs.png", flush=True)         
plt.savefig(f"{plotdir}{sample_name}_hvgs.png") 
plt.close()
adata.raw = adata 
print("Slicing data", flush=True) 
bdata = adata[:, adata.var.highly_variable] 
print("Regressing variables", flush=True) 
bdata = bdata.copy() 
sc.pp.regress_out(bdata, ['total_counts', 'pct_counts_mt']) 
print("Scaling data", flush=True) 
sc.pp.scale(bdata, max_value = 10) 

print(f"Reducing dimensions on {sample_name}", flush=True)
sc.tl.pca(bdata, svd_solver = 'arpack')
sc.pl.pca(bdata, color = ['PHOX2B','PRRX1'], show = False)
print(f"Save fig to {plotdir}{sample_name}_pca_phox2b_prrx1.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_pca_phox2b_prrx1.png")
plt.close()
sc.pl.pca_variance_ratio(bdata, log = True, show = False)
print(f"Save fig to {plotdir}{sample_name}_pca_variance.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_pca_variance.png")
plt.close() 
sc.pp.neighbors(bdata, n_neighbors = 10, n_pcs = 40)
sc.tl.umap(bdata)
sc.pl.umap(bdata, color = ['PHOX2B', 'PRRX1'], show = False)
print(f"Save fig to {plotdir}{sample_name}_umap_phox2b_prrx1.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_umap_phox2b_prrx1.png")
plt.close()
print(f"Clustering {sample_name} data", flush=True)
sc.tl.leiden(bdata, resolution = 0.2)
sc.pl.umap(bdata, color = 'leiden', show = False)
print(f"Save fig to {plotdir}{sample_name}_umap_clusters.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_umap_clusters.png")
plt.close()
bdata.write_h5ad(
    f"{datadir}{sample_name}_ut_scanpy.h5ad",
    compression = hdf5plugin.FILTERS["zstd"]
)
print(f"Saved clustered {datadir}{sample_name} data", flush=True) 

sc.tl.rank_genes_groups(bdata, 'leiden', method = 'wilcoxon')
sc.pl.rank_genes_groups(bdata, n_genes = 25, sharey = False, show = False)
print(f"Save fig to {plotdir}{sample_name}_cluster_markers.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_cluster_markers.png")
plt.close()
sc.tl.paga(bdata, groups = 'leiden')
sc.pl.paga(bdata, show = False)
print(f"Save fig to {plotdir}{sample_name}_paga.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_paga.png")
plt.close()
sc.pl.paga(bdata, color=['PHOX2B', 'PRRX1'], show = False)
print(f"Save fig to {plotdir}{sample_name}_paga_phox2b_prrx1.png", flush=True) 
plt.savefig(f"{plotdir}{sample_name}_paga_phox2b_prrx1.png")
plt.close()

print("Preparing data for MuTrans", flush=True)
sc.tl.tsne(bdata, n_pcs = 30)
sc.pp.neighbors(bdata, metric = 'cosine', n_neighbors = 60, use_rep = 'X')
bdata.write_h5ad(
    f"{datadir}{sample_name}_ut_scanpy_nn.h5ad",
    compression = hdf5plugin.FILTERS["zstd"]
)
print(f"Saved input {datadir}{sample_name} data for MuTrans", flush=True)

print("FINISHED DATA PREP SCRIPT") 