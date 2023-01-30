# papplain

**`papplain`** is a simple R-based pipeline for single cell RNA-seq processing with a focus on integration. `papplain` follows the practices outlined in the [OSCA book](https://bioconductor.org/books/release/OSCA/), with some additional options for integration/batch effect correction methods. 

As a one-stop solution, this package tends to make choices in lieu of the users, with the proviso that these choices follow either defaults or sensible implementations. However, this means that a certain degree of freedom is removed from the end user. This philosophy posits that users who desire total control on the process (or granular specification of parameters) do not need `papplain` and would be more comfortable setting up their own pipelines. Thus `papplain` exists to automate routine analyses the way I usually do them, and offer "quick and dirty" access for exploratory data analysis.

*Where does the name "papplain" come from?* it's an inside joke with some friends and ex colleagues.

# Install

For the time being, clone the repo and install:

```{r}
devtools::install("path/to/cloned/git")
```

The package will require a number of `BioConductor` dependencies. You can install them as follows:

```{r}
BiocManager::install(c("scran", "scuttle", "bluster", "scater", "batchelor", "DropletUtils", "AUCell", "harmony"))
```

# Usage

The pipelines in `papplain` assume that you are working with the output of CellRanger (or something similar) and you imported it into a `SingleCellExperiment` object (hereafter SCE) using `DropletUtils::read10Xcounts()`. This is relevant for gene identifiers, since the `rowData` slot of the SCE will have a "Symbol" and a "ID" column. 

For demo purposes we can use a publicly available dataset, XXXXXX et al. 20XX, which we retrieve using the `scRNAseq` package:

```{r}
sce <- scRNAseq::SegerstolpePancreasData()
colnames(rowData(sce)) = c("Symbol", "ID")
```
Assuming that you have formed a SCE object containing the "individual" column that identifies different batches, you can run an integration pipeline as follows:

```{r}
sce <- papplain(sce, name = "myproject", batch = "individual", verbose = TRUE, save_plots = FALSE)
```

The `name` argument defines the name of the folder that will be created to store files and plots. We set `verbose = TRUE` to print the progress of the pipeline.

The `papplain()` function is a wrapper around a few modules or sub-pipelines that have different degrees of customization. 

The scheme is:

```
papplain
    ├── Quality Control [QC]
    |    ├── run emptyDrops (optional) [QC/EMPTY]
    |    ├── score mito/ribo/malat1 subsets (optional)
    |    ├── filter out (optional)
    |    └── doublet finding (optional) [QC/DBL]
    ├── Normalization and dimensionality reduction [NOR]
    |    ├── pre-clustering
    |    ├── computing pooled factors
    |    ├── log-normalization (simple or multi-batch)
    |    ├── HVG finding (simple or multi-batch)
    |    ├── PCA
    |    └── UMAP
    ├── Integration [INT] (optional) - choose one method
    |    ├── fastMNN [INT/MNN]
    |    |    ├── integration
    |    |    └── UMAP
    |    ├── Seurat [INT/SEURAT]
    |    |    ├── conversion to Seurat
    |    |    ├── normalization and HVG finding
    |    |    ├── find integration anchors
    |    |    ├── integrate data
    |    |    ├── scale data
    |    |    ├── PCA
    |    |    └── UMAP
    |    ├── LIGER [INT/LIGER]
    |    |    ├── conversion to LIGER
    |    |    ├── normalization
    |    |    ├── HVG finding
    |    |    ├── scale data
    |    |    ├── NMF
    |    |    ├── quantile normalization
    |    |    └── UMAP
    |    ├──▷ Harmony [INT/HARMONY]
    |    |    ├── Harmony matrix (on PCA)
    |    |    └── UMAP
    |    ├── Regression [INT/regression]
    |    |    ├── regression on logcounts
    |    |    ├── PCA
    |    |    └── UMAP
    ├── Cell type annotation (optional/in development)
              ├── rank genes
              ├── AUCell
              └── best score assignment
```

Most of the choices can be made around the integration method. `papplain` has implemented 5 methods: `fastMNN` and `regressBatches` from the `batchelor` package, `Harmony`, the CCA-based `Seurat` method, and non-negative matrix factorization (NMF) from `LIGER` through the `rliger` package. 

`LIGER` and `Seurat` integration methods require an intermediate step where package-specific objects are created and some pre-processing steps are repeated again according to the best practices published by the authors of those packages. 

Each step of the pipeline can be called independently on the object:

```{r}
sce <- doQC(sce, name = "segerstolpe", batch = "individual", save_plots = FALSE)
sce <- doNormAndReduce(sce, name = "segerstolpe", batch = "individual")
sce <- integrateSCE(sce, batch = "individual", method = "Seurat")
```

## Parallelization

Since `papplain` is mostly based on R/Bioconductor packages, it offers the `BiocParallel` parallelization backend through its `parallel_param` argument. In some cases, e.g. the Seurat integration, another type of backend has to be set up separately outside of the function call. 

BiocParallel parallelization is implemented where possible, i.e. all the steps of the pipeline where it is sensible to use them (PCA, clustering, integration...). 

```{r}
sce <- integrateSCE(sce, batch = "individual", method = "fastMNN", parallel_param = MulticoreParam(workers = 4))
```

# Clustering

`papplain` also offers a wrapper around clustering functions. For now, only SNN-based Louvain and Leiden clustering are implemented. The `makeGraphsAndClusters()` function allows users to do parameter sweeps along either the number of neighbors or the resolution of the clustering. 

In this example we sweep along the value of the `resolution` parameter for a Louvain clustering. For the SNN graph constructions, edges are weighted according to the jaccard index of their shared neighbors, mimicking the `Seurat` graph construction and clustering procedure:

```{r}
sce <- makeGraphsAndClusters(sce, k = c(0.1, 0.25, 0.5, 0.75, 1),
                            sweep_on = "clustering", method = "louvain", 
                            weighting_scheme = "jaccard", prefix = "SNN_",
                            verbose = TRUE)
```

If another integration method has been run on the same object (e.g. Seurat integration), then the clustering can be performed on that integrated space by specifying the `space` argument (in this case, `reducedDim(sce, "PCA_Seurat")`):

```
sce <- makeGraphsAndClusters(sce, k = c(0.1, 0.25, 0.5, 0.75, 1), 
                            space = reducedDim(sce, "PCA_Seurat")
                            sweep_on = "clustering", method = "louvain", 
                            weighting_scheme = "jaccard", prefix = "Seurat_SNN_",
                            verbose = TRUE)
```

The default value for `space` is `NULL` and will use the `"PCA"` slot from the `reducedDim()` accessor.  

## Plotting clustering results

If using `clustree`, the clustering tree can be visualized by using the same prefix defined in `makeGraphsAndClusters()`:

```{r}
library(clustree)
clustree(sce, prefix = "SNN_")
```

The default arguments to the clustering wrapper include the generation of modularity and approximate silhouette scores for every clustering round. These will be stored in the `metadata` of the SCE, named according to the prefix, the resolution, and the `"modularity_"` and `"silhouette_"` prefixes. 
Silhouette and modularity can be visualized by using the dedicated functions:

```{r}
plotSilhouette(sce, "SNN_0.82")

plotModularity(sce, "SNN_0.82")
```

## Metaclusters
Aditionally, *metaclusters* can also be identified. A metacluster is a cluster of clusters obtained by different clustering methods. Clusters across methods are linked acording to how many cells they share, and these links become edges of a graph. Then, Louvain clustering is run on the graph and the communities that are identified are metaclusters. These metaclusters show the relationship between clustering methods. Moreover, they can be used to understand cluster stability along different parameters and/or integration methods. A cell can belong to different clusters according to the clustering method (i.e. to the resolution or to the space that was used). If a cell belongs to clusters that are consistently included in a metacluster, then that cell belongs to a "stable" cluster. If instead the cell belongs to clusters that have different metacluster assignments, then it's in an "unstable" position, meaning it may be clustered differently according to integration methods and/or resolutions.

The `metaClusters()` function takes a `clusters` argument, which is a vector of column names from the colData of the SCE where clustering results are stored. In this example it is easy to isolate by using `grep()` and searching for the prefix "SNN_".

```{r}
clusterlabels <- colnames(colData(sce))[grep("SNN_", colnames(colData(sce)))]
sce <- metaCluster(sce, clusters = clusterlabels)
```

Every cell will belong to a series of clusters, which in turn belong to a metacluster. For every cell, we count how many times they are assigned to a particular metacluster, and the maximum metacluster is assigned, together with a "metacluster score" (i.e. the frequency of assignment to the maximum metacluster) and whether this score is above or below a certain threshold (0.5 by default). These columns are saved in the colData slot of the SCE.

