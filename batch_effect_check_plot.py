import argparse
import os
from anndata import read_h5ad
import scanpy as sc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--figdir")
args = parser.parse_args()

os.makedirs(args.figdir, exist_ok=True)

adata = read_h5ad(args.input)

sc.pp.pca(adata)

# Compute neighborhood graph
sc.pp.neighbors(adata)

# Run UMAP
sc.tl.umap(adata)

# Plot UMAP
color_list = ["tab:pink", "tab:cyan", "tab:blue", "tab:orange"]
fig = sc.pl.umap(adata, color="unique_core_id", palette=color_list, return_fig=True)
fig.tight_layout()

fig.savefig(f"{args.figdir}/batch_umap.png")
