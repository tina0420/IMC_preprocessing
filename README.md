# Preprocessing for Image Mass Cytometry dataset
## Updated: November 29, 2025. By Tina (Yi-Ting) Hsu 
### Background and Rationale
**IMC_preprocess** is a workflow for preprocessing [Imaging Mass Cytometry™](https://www.fluidigm.com/applications/imaging-mass-cytometry) (IMC) dataset. IMC is a spatial modality with single-cell proteiomics expression data as well as spatial context represented by x-y coordinates of all the cells in images. In the analysis work of IMC, we usually receive the dataset in batches due to the generation time (it is commonly 5-6 hours per core). Hence, it costs more and more time on repeating preprocessing steps since we need to change all the parameters everytime for the new dataset. To adress this, this pipeline helps execute all the preprocessing works before cell type annotation in an automatic way which also establishes standardization in the projet. This workflow includes the typical preprocessing steps in the field, such as quality control and clustering, and the specific procedure of the workflow is documented below. Other single cell expression data with changes to the input data format is available as well for the pipeline which is described below in detail. 

The underlying approach to this workflow has 3 modules:
- Batch Effect Checking: Show the batch effect situation by [UMAP plot](https://umap-learn.readthedocs.io/en/latest/). 
- Leiden Clustering: Run Leiden Clustering to make all the cells in clusters/categories.
- Single-cell Protein Expression across Clusters: Produce a protein expression heatmap of each cell across clusters.

<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/DAG.png?raw=true" width="300" height="300">

Package Dependencies:
*You should have Python 3 and Miniconda installed for Linux*
  - nextflow
  - matplotlib = 3.3.3
  - numpy = 1.19.4
  - pandas = 1.1.4
  - seaborn = 0.11.0
  - pip = 20.2.4

### Usage
The execution of the package *requires* a **Linux OS**. If you are using Windows, you can also run it by **wsl**.

1. Clone this file in a directory you prefer. This will download all the files needed for running the pipeline. The output will be stored in the folder created as "output".
```
git clone https://github.com/
```
2. Create an environment with all the dependencies. 
```
conda env create --name nf_imc_preprocessing --file env
```
3. Activate the environment.
```
conda activate nf_imc_preprocessing
```
4. Now run IMC_preprocess by executing the following code. This will create two files in the **output** folder: `batch_umap.png`, `new_adata.h5ad`, `cluster_result.png`.
```
nextflow run workflow.nf
```
5. Generate a DAG from the workflow:
```
snakemake --dag results/heatmap.png | dot -Tsvg > dag.svg
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

- **new_adata.h5ad**: an anndata object with a new column, leiden in obs, which represents the labels/clusters from Leiden Clustering.

------------------------------------------------------
<img src="https://github.com/tina0420/IMC_preprocessing/blob/main/output/cluster_result.png?raw=true" width="300" height="300">

- **cluster_result.png**: This graph shows the clustering results by the marker expressions ofr each cell across clusters. The objective of the plot is to help finish cell type annotation.For each cluster, we can give the cell type based on the specific higher expressions in the corresponding protein markers/ For example, if most of the cells in Cluster A highly express in CD20 that we may annotate the cluster as B Cell cluster.




[+] Building 2.1s (5/6)                                    docker:desktop-linux
 => [internal] load build definition from Dockerfile                       0.0s
 => => transferring dockerfile: 680B                                       0.0s
 => [internal] load metadata for docker.io/library/python:3.10             0.4s
 => [internal] load .dockerignore                                          0.0s
 => => transferring context: 2B                                            0.0s
 => CACHED [1/3] FROM docker.io/library/python:3.10@sha256:20ca17b2908b02  0.0s
 => => resolve docker.io/library/python:3.10@sha256:20ca17b2908b0202fc975  0.0s
 => ERROR [2/3] RUN apt-get update && apt-get install -y --no-install-rec  1.6s
------                                                                          
 > [2/3] RUN apt-get update && apt-get install -y --no-install-recommends     build-essential     python3-dev     libhdf5-dev     libatlas-dev     gfortran     git     wget     curl     ca-certificates     pkg-config     && apt-get clean     && rm -rf /var/lib/apt/lists/*:                                               
0.126 Hit:1 http://deb.debian.org/debian trixie InRelease                       
0.127 Get:2 http://deb.debian.org/debian trixie-updates InRelease [47.3 kB]
0.137 Get:3 http://deb.debian.org/debian-security trixie-security InRelease [43.4 kB]
0.147 Get:4 http://deb.debian.org/debian trixie/main arm64 Packages [9607 kB]
0.306 Get:5 http://deb.debian.org/debian trixie-updates/main arm64 Packages [5404 B]
0.313 Get:6 http://deb.debian.org/debian-security trixie-security/main arm64 Packages [76.3 kB]
0.965 Fetched 9779 kB in 1s (11.3 MB/s)
0.965 Reading package lists...
1.234 Reading package lists...
1.498 Building dependency tree...
1.569 Reading state information...
1.578 Package libatlas-dev is not available, but is referred to by another package.
1.578 This may mean that the package is missing, has been obsoleted, or
1.578 is only available from another source
1.578 
1.579 E: Package 'libatlas-dev' has no installation candidate
------
Dockerfile:6
--------------------
   5 |     # Install system dependencies
   6 | >>> RUN apt-get update && apt-get install -y --no-install-recommends \
   7 | >>>     build-essential \
   8 | >>>     python3-dev \
   9 | >>>     libhdf5-dev \
  10 | >>>     libatlas-dev \
  11 | >>>     gfortran \ 
  12 | >>>     git \
  13 | >>>     wget \ 
  14 | >>>     curl \
  15 | >>>     ca-certificates \
  16 | >>>     pkg-config \
  17 | >>>     && apt-get clean \
  18 | >>>     && rm -rf /var/lib/apt/lists/*
  19 |     
--------------------
ERROR: failed to build: failed to solve: process "/bin/sh -c apt-get update && apt-get install -y --no-install-recommends     build-essential     python3-dev     libhdf5-dev     libatlas-dev     gfortran     git     wget     curl     ca-certificates     pkg-config     && apt-get clean     && rm -rf /var/lib/apt/lists/*" did not complete successfully: exit code: 100

