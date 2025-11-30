import argparse
import os
from anndata import read_h5ad
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--output")
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

adata = read_h5ad(args.input)

# Run leiden clustering
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.leiden(adata)

adata.write(f"{args.output}/new_adata.h5ad")


