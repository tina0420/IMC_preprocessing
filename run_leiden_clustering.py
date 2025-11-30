import argparse
import os
from anndata import read_h5ad
import scanpy as sc
from sklearn.mixture import BayesianGaussianMixture

parser = argparse.ArgumentParser()
parser.add_argument("--input")
parser.add_argument("--output")
args = parser.parse_args()

adata = read_h5ad(args.input)

# Run leiden clustering
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, use_rep='X_pca')

# Run PCA
sc.tl.pca(adata, svd_solver='arpack')

# Extract PCA features for clustering
X_pca = adata.obsm['X_pca']

# Run Bayesian Gaussian Mixture Model
bgmm = BayesianGaussianMixture(
    n_components=10,
    covariance_type='full', 
    weight_concentration_prior_type='dirichlet_process', 
    max_iter=1000, 
    random_state=0
)
adata.obs['bgmm'] = bgmm.fit_predict(X_pca).astype(str)

adata.write(f"{args.output}/new_adata.h5ad")


