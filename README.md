# Preprocessing for Imaging Mass Cytometry dataset
## Updated: November 29, 2025. By Tina (Yi-Ting) Hsu 
### Background and Rationale
**IMC_preprocess** is a workflow for preprocessing [Imaging Mass Cytometry™](https://www.fluidigm.com/applications/imaging-mass-cytometry) (IMC) dataset. IMC is a spatial modality with single-cell proteiomics expression data as well as spatial context represented by x-y coordinates of all the cells in images. In the analysis work of IMC, we usually receive the dataset in batches due to the generation time (it is commonly 5-6 hours per core). Hence, it costs more and more time on repeating preprocessing steps since we need to change all the parameters everytime for the new dataset. To address this, this pipeline helps execute all the preprocessing works before cell type annotation in an automatic way which also establishes standardization in the projet. This workflow includes the typical preprocessing steps in the field, such as quality control and clustering, and the specific procedure of the workflow is documented below. Other single cell expression data with changes to the input data format is available as well for the pipeline which is described below in detail. 

The underlying approach to this workflow has 3 modules:
- Batch Effect Checking: Show the batch effect situation by [UMAP plot](https://umap-learn.readthedocs.io/en/latest/). 
- Bayesian Gaussian Mixture Model Clustering: Run Bayesian GMM Clustering to make all the cells in clusters/categories.
- Single-Cell Protein Expression across Clusters: Produce a protein expression heatmap of each cell across clusters.

<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/flowchart.png?raw=true" width="300" height="300">

Package Dependencies:
*You should have Python 3 and Miniconda installed for Linux*
  - nextflow
  - matplotlib = 3.10.1
  - numpy = 2.1.3
  - pandas = 2.2.3
  - seaborn = 0.13.2
  - scikit-learn = 1.5.2
  - anndata = 0.11.4
  - scanpy = 1.11.1

### Usage
The execution of the package *requires* a **Linux OS**. If you are using Windows, you can also run it by **wsl**.

1. Clone this file in a directory you prefer. This will download all the files needed for running the pipeline. The output will be stored in the folder created as "output".
```
git clone https://github.com/tina0420/IMC_preprocessing.git
```
2. Create an environment with all the dependencies. 
```
conda env create --name imc_preprocessing --file environment.yaml
```
3. Activate the environment.
```
conda activate imc_preprocessing
```
4. Now run IMC_preprocess by executing the following code. This will create two files in the **figures** folder: `batch_umap.png`, `new_adata.h5ad`, `cluster_result.png`.
```
nextflow run workflow.nf
```
5. Alternatively, generate a DAG along with executing the workflow:
```
nextflow run workflow.nf -with-dag flowchart.png
```

### Input
IMC_preprocessing takes in one h5ad file: an anndata-structure file with all quantified spatial data. All the files should be placed in the folder `data/`. 

Inside `data/`, the sample dataset we use here is a segmented and quantified IMC dataset with 4 breast cancer patients represented. All protein expressions and x-y coordinates of per cell should be stored in the h5ad.

H5AD should have the following properties:
*If you are using the spatial data from other technologies, please make sure the data type/format is h5ad and the following properties stored in the correct layers in an anndata object. *

- anndata.X: expression values of all the markers
- anndata.obsm["spatial"]: x-y coordinates of all the cells across iamges.
- anndata.obs["unique_core_id"]: the column in obs with the unique key to the cores or images.

### Output
IMC_preprocessing produces two types of visualizations in png files for interpreting data characteristics in preprocessing step: UMAP plot, expression heatmap plot, and one h5ad file with the clustering results. All the files should be stored in the folder `output/`.

Showcase for the output from IMC_preprocessing 
<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/output/batch_umap.png?raw=true" width="300" height="300">

- **batch_umap.png**: This is a visualization for high dimensional data based on Uniform Manifold Approximation and Projection. Axes do not refer to spatial coordinates. A UMAP plot is included to show that how different are the images, cores or batches. This is a common step in quality control that to evaluate if batch correction is needed or not. If there's no series batch effect across batches, then we can run the clustering with the whole adata for annotation.

------------------------------------------------------

- **new_adata.h5ad**: an anndata object with a new column, bgmm in obs, which represents the labels/clusters from BGMM Clustering.

------------------------------------------------------
<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/output/cluster_result.png?raw=true" width="300" height="300">

- **cluster_result.png**: This graph shows the clustering results by the marker expressions for each cell across clusters. The objective of the plot is to help finish cell type annotation.For each cluster, we can give the cell type based on the specific higher expressions in the corresponding protein markers/ For example, if most of the cells in Cluster A highly express in CD20 that we may annotate the cluster as B Cell cluster.

