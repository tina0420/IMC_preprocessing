import argparse
import os
from anndata import read_h5ad
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--figdir")
args = parser.parse_args()

adata = read_h5ad(args.input)

labels = adata.obs["bgmm"].values

df = pd.DataFrame(adata.X, columns=adata.var.index).assign(label=labels)
df = df.sort_values(by="label")

color_list = [
    "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#800000", "#008000",
    "#000080", "#808000", "#800080", "#008080", "#C00000", "#00C000", "#0000C0", "#C0C000",
    "#C000C0", "#00C0C0", "#FF8000", "#FF0080", "#80FF00", "#00FF80", "#0080FF", "#8000FF",
    "#FF8080", "#80FF80", "#8080FF", "#FF80FF", "#80FFFF", "#FFFF80", "#404040", "#808080",
    "#C0C0C0", "#2020FF", "#20FF20", "#FF2020", "#20FFFF", "#FF20FF", "#FFFF20"
]

if len(pd.unique(df["label"])) > len(color_list):
    color_list += color_list

color_mapping = dict(zip(sorted(list(pd.unique(df["label"]))), range(len(pd.unique(df["label"])))))
row_colors = pd.DataFrame({'Cluster': [color_list[color_mapping[i]] for i in df["label"]]})

df = df.drop(columns=["label"]).reset_index(drop=True)

g = sns.clustermap(df, row_colors=row_colors, row_cluster=False, col_cluster=False, figsize=(8, 8), cmap=sns.diverging_palette(240, 10, n=9), cbar_pos=None)

ax = g.ax_heatmap
plt.xticks(fontsize=10, rotation=90, ha="right")
ax.set_yticklabels([])
ax.tick_params(right=False, length=1, width=0.5)
g.ax_row_colors.set_xticklabels(["Cluster"], fontsize=8, rotation=90, ha="right")
g.ax_row_colors.tick_params(length=1, width=0.5)

plt.title("Cell Expression", weight='bold', fontsize=15)
plt.tight_layout()

plt.savefig(f"{args.figdir}/cluster_result.png")
plt.clf()
