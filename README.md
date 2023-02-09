# papplain

**`papplain`** is a simple R-based pipeline for single cell RNA-seq processing with a focus on integration.

`papplain` follows the practices outlined in the [OSCA book](https://bioconductor.org/books/release/OSCA/), with some additional options for integration/batch effect correction methods.

As a one-stop solution, this package tends to make choices for the users, with the caveat that these choices follow either defaults or sensible implementations. However, this means that a certain degree of freedom is removed from the end user. This assumes that users who desire total control on the process (or granular specification of parameters) do not need `papplain` and would be more comfortable setting up their own analysis pipelines.

`papplain` exists to automate and share routine analyses the way I usually do them, and offer "quick and dirty" access for exploratory data analysis.

*Where does the name "papplain" come from?* it's an inside joke with some friends and ex colleagues.

# Install

For the time being, clone the repo and install:

```{r}
devtools::install("path/to/cloned/git")
```

The package will require a number of `BioConductor` and `GitHub` dependencies. You can install them as follows:

```{r}
BiocManager::install(c("scran", "scuttle", "bluster", "scater", "batchelor", "DropletUtils", "AUCell", "harmony", "GSVA",  "gdagstn/alamak"))
```

# Usage

The pipelines in `papplain` assume that you are working with the output of CellRanger (or something similar) and you imported it into a `SingleCellExperiment` object (hereafter SCE) using `DropletUtils::read10Xcounts()`. This is relevant for gene identifiers, since the `rowData` slot of the SCE will have a "Symbol" and a "ID" column.

For demo purposes we can use a publicly available dataset, Segerstolpe et al. 2016 [ref](https://pubmed.ncbi.nlm.nih.gov/27667667/), which we retrieve using the `scRNAseq` package:

```{r}
# BiocManager::install("scRNAseq") 
sce <- scRNAseq::SegerstolpePancreasData()
colnames(rowData(sce)) = c("Symbol", "ID")
```

Assuming that you have formed a SCE object containing the "individual" column that identifies different batches, you can run an integration pipeline as follows:

```{r}
sce <- papplain(sce, name = "myproject", batch = "individual", 
                integration_method = "Harmony",
                verbose = TRUE, save_plots = FALSE)
```

The `name` argument defines the name of the folder that will be created to store files and plots. We set `verbose = TRUE` to print the progress of the pipeline.

The `papplain()` function is a wrapper around a few modules or sub-pipelines that have different degrees of customization.

The scheme is:

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
        |    ├── Harmony [INT/HARMONY]
        |    |    ├── Harmony matrix (on PCA)
        |    |    └── UMAP
        |    └── Regression [INT/regression]
        |         ├── regression on logcounts
        |         ├── PCA
        |         └── UMAP
        └── Cell type annotation [ANNO] (optional) - choose one method
                  ├── Seurat AddModuleScore
                  ├── ssGSEA
                  ├── UCell
                  └── AUCell
                 

Most of the choices can be made around the integration method. `papplain` has implemented 5 methods: `fastMNN` and `regressBatches` from the `batchelor` package, `Harmony`, the CCA-based `Seurat` method, and non-negative matrix factorization (NMF) from `LIGER` through the `rliger` package.

`LIGER` and `Seurat` integration methods require an intermediate step where package-specific objects are created and some pre-processing steps are repeated again according to the best practices published by the authors of those packages.

Each step of the pipeline can be called independently on the object:

```{r}
sce <- doQC(sce, name = "segerstolpe", batch = "individual", integration_method = "Harmony", save_plots = FALSE)
sce <- doNormAndReduce(sce, name = "segerstolpe", batch = "individual")
sce <- integrateSCE(sce, batch = "individual", method = "Seurat")
```

## Plotting

Plotting functions are a work in progress, but for now you can plot a UMAP with a few visualization options.

You can choose point color using the `color_by` argument, and facetting is supported via the `group_by` argument. Additionally you can choose a `shape_by` for symbols, and `label_by` to place labels on the plot. Note that shape, group, and labels need to be categorical (i.e. factor) variables, whereas color can be numeric. The color palette is automatically generated, but it can be set by the user through the `color_palette` argument.

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "individual", group_by = "disease")
```

<img src="https://user-images.githubusercontent.com/21171362/216002545-43210b08-5919-49ce-8688-e42b8bf70a64.png" width="800"/>

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "sum")
```

<img src="https://user-images.githubusercontent.com/21171362/216003504-32242a1e-bf33-4479-b525-c29a6569d64f.png" width="300"/>

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
                             space = reducedDim(sce, "PCA_Harmony")[,1:20],
                             sweep_on = "clustering", method = "louvain", 
                             weighting_scheme = "jaccard", prefix = "SNN_",
                             verbose = TRUE)
```

If another integration method has been run on the same object (e.g. Seurat integration), then the clustering can be performed on that integrated space by specifying the `space` argument (in this case, `reducedDim(sce, "PCA_Seurat")`):

```{r}
sce <- makeGraphsAndClusters(sce, k = c(0.1, 0.25, 0.5, 0.75, 1), 
                             space = reducedDim(sce, "PCA_Seurat")
                             sweep_on = "clustering", method = "louvain", 
                             weighting_scheme = "jaccard", prefix = "Seurat_SNN_",
                             verbose = TRUE)
```

The default value for `space` is `NULL` and will use the `"PCA"` slot from the `reducedDim()` accessor.

## Plotting clustering results

You can visualize clustering results on the UMAP using the `plot_UMAP()` function, adding labels if desired:

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "SNN_0.5", label_by = "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216006173-db0a4ebf-702e-4bac-a25c-e789ee675e03.png" width="400"/>

If using `clustree`, the clustering tree can be visualized by using the same prefix defined in `makeGraphsAndClusters()`:

```{r}
library(clustree)
clustree(sce, prefix = "SNN_")
```

<img src="https://user-images.githubusercontent.com/21171362/216004500-88ec933a-3fc4-486b-bdfa-ba5aae2092ad.png" width="329"/>

The default arguments to the clustering wrapper include the generation of modularity and approximate silhouette scores for every clustering round. These will be stored in the `metadata` of the SCE, named according to the prefix, the resolution, and the `"modularity_"` and `"silhouette_"` prefixes. Silhouette and modularity can be visualized by using the dedicated functions:

```{r}
plotSilhouette(sce, "SNN_0.5")

plotModularity(sce, "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216005004-41ae56a0-adae-40f3-8154-6ce59ae7a3ed.png" width="400"/> <img src="https://user-images.githubusercontent.com/21171362/216005131-d7a639e1-928c-4dcd-b695-573427a4d14f.png" width="400"/>

## Gene dot plot

You can also use the `plot_dots()` function to plot the popular dot-plot for marker genes.

This function takes in a `SingleCellExperiment` object, together with a vector of genes (matched to the `rownames` of the object), and a grouping variable specified by the `group_by` argument. Additionally, dots can be ordered by hierarchical clustering on either genes, groups, or both (set `cluster_genes` and/or `cluster_groups` to `TRUE`, which is the default). Colors can also be customized via the `color_palette` argument. Finally, the user can choose whether they want genes to be columns (`format = "wide"`, the default) or rows (`format = "tall"`).

```{r}
# Quick and dirty marker calculation
markers = presto::wilcoxauc(sce, group_by = "SNN_0.5")
markerlist = split(markers, markers$group)

for(i in seq_len(length(markerlist))) {
  markerlist[[i]]$deltapct = markerlist[[i]]$pct_in - markerlist[[i]]$pct_out
  markerlist[[i]] = markerlist[[i]][order(markerlist[[i]]$deltapct, decreasing = TRUE),]
}

top5 = lapply(markerlist, function(x) x$feature[1:5])
markergenes =  Reduce(union, top5)

plot_dots(sce, genes = top5, group_by = "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216310875-a06081b7-e9bf-404a-99f3-ef1291555ad3.png" width="600"/>


## Metaclusters

Aditionally, *metaclusters* can also be identified. A metacluster is a cluster of clusters obtained by different clustering methods. Clusters across methods are linked acording to how many cells they share, and these links become edges of a graph. Then, Louvain clustering is run on the graph and the communities that are identified are metaclusters. These metaclusters show the relationship between clustering methods. Moreover, they can be used to understand cluster stability along different parameters and/or integration methods. A cell can belong to different clusters according to the clustering method (i.e. to the resolution or to the space that was used). If a cell belongs to clusters that are consistently included in a metacluster, then that cell belongs to a "stable" cluster. If instead the cell belongs to clusters that have different metacluster assignments, then it's in an "unstable" position, meaning it may be clustered differently according to integration methods and/or resolutions.

The `metaClusters()` function takes a `clusters` argument, which is a vector of column names from the colData of the SCE where clustering results are stored. In this example it is easy to isolate by using `grep()` and searching for the prefix "SNN\_".

```{r}
clusterlabels <- colnames(colData(sce))[grep("SNN_", colnames(colData(sce)))]
sce <- metaCluster(sce, clusters = clusterlabels)
```

<img src="https://user-images.githubusercontent.com/21171362/216005406-369d20fb-5696-4d91-ba0a-b5e450d8db86.png" width="400"/>

Every cell will belong to a series of clusters, which in turn belong to a metacluster. For every cell, we count how many times they are assigned to a particular metacluster, and the maximum metacluster is assigned, together with a "metacluster score" (i.e. the frequency of assignment to the maximum metacluster) and whether this score is above or below a certain threshold (0.5 by default). These columns are saved in the `colData` slot of the SCE.

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "metacluster_score", label_by = "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216006433-2b39bf37-a9f4-49e4-be97-ec17ac690297.png" width="400"/>

## Assigning cell identities

`papplain` implements two methods for automated cell identity assignment, based on the Bioconductor `AUCell` package, the `GSVA` `ssGSEA` implementation, or the `Seurat` `AddModuleScore()` function.

The function requires a `genesets` named list containing genes to be used for scoring every single cell. These can be obtained through other packages, e.g. `msigdbr`. For instance, if we wanted to take all the Muraro et al. signature genes, present in the C8 collection, we would do:

```{r}
library(msigdbr)

type_genes = msigdbr("Homo sapiens", category = "C8")
genesets = lapply(split(type_genes, type_genes$gs_name), function(x) x$gene_symbol)
muraro_genes = genesets[grep("MURARO", names(genesets))]
```

Then, we would use the `assignIdentities()` function from `papplain` to calculate signature scores:

```{r}
sce = assignIdentities(sce, 
                       genesets = muraro_genes, 
                       method = "AUC")
```

Other methods are "Seurat", "UCell", and "ssGSEA". 

This will create a column named "labels_AUC" (or anything else the user determines using the `name` argument) in the `colData(sce)`. Assignments can be plotted:

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "labels_AUC")
```


<img src="https://user-images.githubusercontent.com/21171362/217058249-6bcc821a-22cb-4aa1-88d1-d21e33fc63e7.png" width = "800"/>