<img src="https://user-images.githubusercontent.com/21171362/232594749-6a6beacc-6c59-4e40-ae62-e7026f355742.png" align="right" width="200"/>

# cellula

**`cellula`** is a simple R-based pipeline for single cell RNA-seq
processing with a number of methods for integration and identity
assignment.

`cellula` follows the practices outlined in the [OSCA
book](https://bioconductor.org/books/release/OSCA/)[[1]](#1), with some
additional options for integration/batch effect correction methods,
signature scoring, and downsampling.

As a one-stop solution, this package tends to make choices for the
users, with the caveat that these choices follow either defaults or
sensible implementations. However, this means that a certain degree of
freedom is removed from the end user. This assumes that users who desire
total control on the process (or granular specification of parameters)
do not need `cellula` and would be more comfortable setting up their own
analysis pipelines.

`cellula` exists to automate and share routine analyses the way I
usually do them, and offer "quick and dirty" access for exploratory data
analysis.

`cellula` is very much under **active development**, and any feedback
and contribution are welcome through the
[Issues](https://github.com/gdagstn/cellula/issues) page.

*Where did the name "papplain" come from?* it's an inside joke with some
friends and ex colleagues.

*Why did you change it to "cellula"?* one day I'd like to share this
tool and I need a name that is not too dumb.

# Table of Contents

-   [Install](https://github.com/gdagstn/cellula#install)
-   [Usage](https://github.com/gdagstn/cellula#usage)
    -   [Plotting](https://github.com/gdagstn/cellula#plotting)
    -   [Parallelization](https://github.com/gdagstn/cellula#parallelization)
-   [Clustering](https://github.com/gdagstn/cellula#clustering)
    -   [Plotting clustering
        results](https://github.com/gdagstn/cellula#plotting-clustering-results)
    -   [Gene dot
        plot](https://github.com/gdagstn/cellula#gene-dot-plot)
    -   [Metaclusters](https://github.com/gdagstn/cellula#metaclusters)
-   [Assigning cell
    identities](https://github.com/gdagstn/cellula#assigning-cell-identities)
-   [Downsampling](https://github.com/gdagstn/cellula#downsampling)
-   [Inferring
    trajectories](https://github.com/gdagstn/cellula#inferring-trajectories)
    -   [Metacells](https://github.com/gdagstn/cellula#metacells)

# Install

Use `remotes`, `devtools` or `BiocManager` to install:

```{r}
remotes::install_github("gdagstn/cellula")
devtools::install_github("gdagstn/cellula")
BiocManager::install("gdagstn/cellula")

```

`cellula` is **dependency-heavy**, which is not something I'm proud of,
but makes sense considering this is a wrapper to a series of different
analytical approaches.

The package will require a number of `BioConductor` and `GitHub`
dependencies. You can install them as follows:

```{r}
BiocManager::install(c("scran", "scuttle", "bluster", "scater", "batchelor", "DropletUtils", 
                      "AUCell", "harmony", "GSVA", "gdagstn/oveRlay", "UCell", 
                      "slingshot", "TSCAN", "cole-trapnell-lab/monocle3"))
```

# Usage

The pipelines in `cellula` assume that you are working with the output
of `CellRanger` (or something similar) and you imported it into a
`SingleCellExperiment` object (hereafter SCE) using the `TENxIO`
package. This is relevant for gene identifiers, since the `rowData` slot
of the SCE will have a "Symbol" and a "ID" column.

For demo purposes we can use a publicly available dataset, [Segerstolpe
et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27667667/)[[2]](#2), which
we retrieve using the `scRNAseq` package:

```{r}
# BiocManager::install("scRNAseq") 
sce <- scRNAseq::SegerstolpePancreasData()
colnames(rowData(sce)) = c("Symbol", "ID")
```

Assuming that you have formed a SCE object containing the "individual"
column that identifies different batches, you can run an integration
pipeline as follows:

```{r}
sce <- cellula(sce, name = "myproject", batch = "individual", 
                integration_method = "Harmony",
                verbose = TRUE, save_plots = TRUE)
```

The `name` argument defines the name of the folder that will be created
to store files and plots. We set `verbose = TRUE` to print the progress
of the pipeline. Setting `save_plots = TRUE` will create a few QC plots
in the `name/plots` folder: total UMI, total genes detected, UMI x
genes; optionally % MT, % Ribo and %MALAT1, total UMI x doublet class.
Plots are separated according to whether the cells were discarded or not
in the filtering step.

The `cellula()` function is a wrapper around a few modules or
sub-pipelines that have different degrees of customization. 

There are other independent functions that are not run through `cellula()` 
as they need some user input, e.g. `findTrajectories()` requires the user
to specify the starting cluster (through `makeGraphsAndClusters()`), or 
the cluster labels to use. 

The scheme is:

```         
cellula()
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
              ├── AUCell
              └── Jaitin
              
  makeGraphsAndClusters()
    └── Multi-resolution clustering [CLU]
          ├── sweep on Louvain/Leiden resolution or SNN neighbor numbers
          ├── calculate modularity (optional)
          └── calculate silhouette (optional)
              
  findTrajectories()            
    └── Trajectory estimation [TRAJ]
          ├── slingshot [TRAJ/slingshot]
          |   ├── get lineages
          |   ├── calculate principal curves
          |   ├── embed in 2D (optional)
          |   └── calculate per-lineage DE (optional)
          └── monocle3 [TRAJ/monocle]
              ├── convert to CellDataSet
              ├── learn graph
              └── embed in 2D (optional)
                   ├── populate FR layout (if dr_embed = "FR")
                   └── UMAP on FR layout (if dr_embed = "FR")
             
```

Most of the choices can be made around the integration method. `cellula`
has implemented 5 methods: `fastMNN`[[3]](#3),[[4]](#4) and
`regressBatches` from the `batchelor` package, `Harmony`[[5]](#5), the
CCA-based `Seurat`[[6]](#6) method, and non-negative matrix
factorization (NMF) from `LIGER`[[7]](#7) through the `rliger` package.

`LIGER` and `Seurat` integration methods require an intermediate step
where package-specific objects are created and some pre-processing steps
are repeated again according to the best practices published by the
authors of those packages.

Each step of the pipeline can be called independently on the object:

```{r}
sce <- doQC(sce, name = "segerstolpe", batch = "individual", save_plots = TRUE)
sce <- doNormAndReduce(sce, name = "segerstolpe", batch = "individual")
sce <- integrateSCE(sce, batch = "individual", method = "Seurat")
```

Doublet identification is carried out through the `scDblFinder`
package[[8]](#8) using standard defaults.

## Plotting

There are a few simple plotting functions in `cellula`:

-   `plot_UMAP()` to plot a UMAP (or any other 2D dimensional reduction)
-   `plot_Coldata()` to plot data from the `colData` slot as a boxplot,
    scatterplot or confusion matrix
-   `plot_dots()` to plot a dot plot of gene expression

### `plot_UMAP()`

You can choose point color using the `color_by` argument, and facetting
is supported via the `group_by` argument. Additionally you can choose a
`shape_by` for symbols, and `label_by` to place labels on the plot. Note
that shape, group, and labels need to be categorical (i.e. factor)
variables, whereas color can be numeric. The color palette is
automatically generated, but it can be set by the user through the
`color_palette` argument.

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "individual", group_by = "disease")
```

<img src="https://user-images.githubusercontent.com/21171362/216002545-43210b08-5919-49ce-8688-e42b8bf70a64.png" width="800"/>

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "sum")
```

<img src="https://user-images.githubusercontent.com/21171362/216003504-32242a1e-bf33-4479-b525-c29a6569d64f.png" width="300"/>

### `plot_Coldata()`

Takes as input `x` and `y` as column names from `colData(sce)`, with an
optional `color_by` and `group_by` argument for facetting.

This function returns different plots depending on the class of the 2
`colData` columns selected: - if `y` is a numeric and `x` is categorical
(character or factor), it returns a combined violin-boxplot with one
plot per level of `x`.

```{r}
plot_Coldata(sce, x = "individual", y = "sum") + scale_y_log10()
```

<img src="https://user-images.githubusercontent.com/21171362/218302530-3e934eaa-20ce-43bf-98a0-b4c211d77ed4.png" width="800"/>

Additionally, if the `color_by` argument specifies another column, every
`x` will be divided by levels of `color_by`. With the appropriate use of
the `x`, `color_by` and `group_by` variables once an look at 3 different
groupings of `y` at once.

```{r}
plot_Coldata(sce, x = "individual", y = "sum", color_by = "disease", group_by = "cell type") + scale_y_log10()
```

-   if `y` and `x` are both categorical, it returns a heatmap of the
    confusion matrix where every value is the pairwise Jaccard index
    between sets for any given level pair (this is mostly useful to
    check for differences in clustering/annotations)

-   if `y` and `x` are both numeric, it returns a scatterplot with an
    optional 2D kernel density contour plot overlaid.

```{r}
plot_Coldata(sce, x = "sum", y = "detected") + scale_x_log10()
```

<img src="https://user-images.githubusercontent.com/21171362/218302610-403b2620-f12c-486f-8055-0dec050f1c55.png" width="500"/>

`plot_dots()` - see below.

## Parallelization

Since `cellula` is mostly based on R/Bioconductor packages, it offers
the `BiocParallel` parallelization backend through its `parallel_param`
argument. In some cases, e.g. the Seurat integration, another type of
backend has to be set up separately outside of the function call.

BiocParallel parallelization is implemented where possible, i.e. all the
steps of the pipeline where it is sensible to use them (PCA, clustering,
integration...).

```{r}
sce <- integrateSCE(sce, batch = "individual", method = "fastMNN", parallel_param = MulticoreParam(workers = 4))
```

# Clustering

`cellula` also offers a wrapper around clustering functions. For now,
only SNN-based Louvain[[9]](#9) and Leiden[[10]](#10) clustering are
implemented. The `makeGraphsAndClusters()` function allows users to do
parameter sweeps along either the number of neighbors or the resolution
of the clustering.

In this example we sweep along the value of the `resolution` parameter
for a Louvain clustering. For the SNN graph constructions, edges are
weighted according to the jaccard index of their shared neighbors,
mimicking the `Seurat` graph construction and clustering procedure:

```{r}
sce <- makeGraphsAndClusters(sce, k = c(0.1, 0.25, 0.5, 0.75, 1),
                             space = reducedDim(sce, "PCA_Harmony")[,1:20],
                             sweep_on = "clustering", method = "louvain", 
                             weighting_scheme = "jaccard", prefix = "SNN_",
                             verbose = TRUE)
```

If another integration method has been run on the same object (e.g.
Seurat integration), then the clustering can be performed on that
integrated space by specifying the `space` argument (in this case,
`reducedDim(sce, "PCA_Seurat")`):

```{r}
sce <- makeGraphsAndClusters(sce, k = c(0.1, 0.25, 0.5, 0.75, 1), 
                             space = reducedDim(sce, "PCA_Seurat"),
                             sweep_on = "clustering", method = "louvain", 
                             weighting_scheme = "jaccard", prefix = "Seurat_SNN_",
                             verbose = TRUE)
```

The default value for `space` is `NULL` and will use the `"PCA"` slot
from the `reducedDim()` accessor.

## Plotting clustering results

You can visualize clustering results on the UMAP using the `plot_UMAP()`
function, adding labels if desired:

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "SNN_0.5", label_by = "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216006173-db0a4ebf-702e-4bac-a25c-e789ee675e03.png" width="400"/>

If using the `clustree` [package]()[[11]](#11), the clustering tree can
be visualized by using the same prefix defined in
`makeGraphsAndClusters()`:

```{r}
library(clustree)
clustree(sce, prefix = "SNN_")
```

<img src="https://user-images.githubusercontent.com/21171362/216004500-88ec933a-3fc4-486b-bdfa-ba5aae2092ad.png" width="329"/>

The default arguments to the clustering wrapper include the generation
of modularity and approximate silhouette scores for every clustering
round. These will be stored in the `metadata` of the SCE, named
according to the prefix, the resolution, and the `"modularity_"` and
`"silhouette_"` prefixes. Silhouette and modularity can be visualized by
using the dedicated functions:

```{r}
plotSilhouette(sce, "SNN_0.5")

plotModularity(sce, "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216005004-41ae56a0-adae-40f3-8154-6ce59ae7a3ed.png" width="400"/>
<img src="https://user-images.githubusercontent.com/21171362/216005131-d7a639e1-928c-4dcd-b695-573427a4d14f.png" width="400"/>

## Gene dot plot

You can also use the `plot_dots()` function to plot the popular dot-plot
for marker genes.

This function takes in a `SingleCellExperiment` object, together with a
vector of genes (matched to the `rownames` of the object), and a
grouping variable specified by the `group_by` argument. Additionally,
dots can be ordered by hierarchical clustering on either genes, groups,
or both (set `cluster_genes` and/or `cluster_groups` to `TRUE`, which is
the default). Colors can also be customized via the `color_palette`
argument. Finally, the user can choose whether they want genes to be
columns (`format = "wide"`, the default) or rows (`format = "tall"`).

```{r}
# Install {presto}
remotes::install_github("immunogenomics/presto")

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

Aditionally, *metaclusters* can also be identified. A metacluster is a
cluster of clusters obtained by different clustering methods. Clusters
across methods are linked acording to how many cells they share, and
these links become edges of a graph. Then, Louvain clustering is run on
the graph and the communities that are identified are metaclusters.
These metaclusters show the relationship between clustering methods.
Moreover, they can be used to understand cluster stability along
different parameters and/or integration methods. A cell can belong to
different clusters according to the clustering method (i.e. to the
resolution or to the space that was used). If a cell belongs to clusters
that are consistently included in a metacluster, then that cell belongs
to a "stable" cluster. If instead the cell belongs to clusters that have
different metacluster assignments, then it's in an "unstable" position,
meaning it may be clustered differently according to integration methods
and/or resolutions.

The `metaClusters()` function takes a `clusters` argument, which is a
vector of column names from the colData of the SCE where clustering
results are stored. In this example it is easy to isolate by using
`grep()` and searching for the prefix "SNN\_".

```{r}
clusterlabels <- colnames(colData(sce))[grep("SNN_", colnames(colData(sce)))]
sce <- metaCluster(sce, clusters = clusterlabels)
```

<img src="https://user-images.githubusercontent.com/21171362/216005406-369d20fb-5696-4d91-ba0a-b5e450d8db86.png" width="400"/>

Every cell will belong to a series of clusters, which in turn belong to
a metacluster. For every cell, we count how many times they are assigned
to a particular metacluster, and the maximum metacluster is assigned,
together with a "metacluster score" (i.e. the frequency of assignment to
the maximum metacluster) and whether this score is above or below a
certain threshold (0.5 by default). These columns are saved in the
`colData` slot of the SCE.

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "metacluster_score", label_by = "SNN_0.5")
```

<img src="https://user-images.githubusercontent.com/21171362/216006433-2b39bf37-a9f4-49e4-be97-ec17ac690297.png" width="400"/>

# Assigning cell identities

`cellula` implements 5 methods for automated cell identity assignment, 4 of which
are signature-based, and one reference-based. The methods are based on the 
Bioconductor `AUCell` [package]()[[12]](#12), the `GSVA` `ssGSEA` 
implementation[[13]](#13), the `Seurat` `AddModuleScore()` function, the `UCell`
method[[14]](#14) and the `Jaitin` method[[15]](#15).

For the first 4 methods (`"AUCell"`, `"ssGSEA"`, `"AddModuleScore"` and `"UCell"`) the function 
requires user-defined `genesets`, i.e. a named list containing genes to be used 
for scoring every single cell. These can be obtained through other packages, 
e.g. `msigdbr`. For instance, if we wanted to take all the Muraro et al.[[16]](#16). 
signature genes, present in the C8 collection, we would do:

```{r}
library(msigdbr)

type_genes = msigdbr("Homo sapiens", category = "C8")

genesets = lapply(split(type_genes, type_genes$gs_name), function(x) x$gene_symbol)

muraro_genes = genesets[grep("MURARO", names(genesets))]
```

Then, we would use the `assignIdentities()` function from `cellula` to
calculate signature scores:

```{r}
sce = assignIdentities(sce, 
                       genesets = muraro_genes, 
                       method = "AUC")
```

Other signature-based methods are `"Seurat"`, `"UCell"`, and `"ssGSEA"`.

A reference-based method, `"Jaitin"`, is available. This method uses a matrix of
gene expression as a reference, and then calculates the posterior probability 
for each cell in the `sce` object that its transcriptome matches any of the 
reference ones. The reference with the highest probability is selected as the 
best label. Importantly, for the "`Jaitin`" method it is possible to choose the
`assay` to be used. If the user supplies a matrix of log-normalized counts as a
reference, the `assay` argument should point to a similarly normalized data,
e.g. `"logcounts"`. Finally, `"Jaitin"` is very slow but can be easily parallelized
by supplying a `BPPARAM` object.

`assignIdentities()` will create a column named "labels_AUC" (or anything else the user
determines using the `name` argument) in the `colData(sce)`. Assignments
can be plotted:

```{r}
plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "labels_AUC")
```

<img src="https://user-images.githubusercontent.com/21171362/217058249-6bcc821a-22cb-4aa1-88d1-d21e33fc63e7.png" width="800"/>

You can now see why it can be useful to plot a confusion matrix as a
heatmap:

```{r}
plot_Coldata(sce, x = "SNN_0.5", y = "labels_AUC")
```

<img src="https://user-images.githubusercontent.com/21171362/218304923-ef02142a-00dc-4065-a4d6-a0d509e88d3a.png" width="800"/>

You can also use single signatures as an input, which will result in
adding the score to the `colData` slot of the SCE directly, rather than
an assignment:

```{r}
sce <- assignIdentities(sce, 
                        genesets = muraro_genes$BETA_CELL, 
                        method = "AUC", 
                        name = "Beta_Cell_signature")

plot_UMAP(sce, umap_slot = "UMAP_Harmony", color_by = "Beta_Cell_signature")
```

<img src="https://user-images.githubusercontent.com/21171362/218026010-9bc1f806-e3c3-408d-ad49-0ca5cab1bd09.png" width="500"/>

The `"UCell"` method works well when you have small signatures (e.g.
even 2/3 genes). It allows you to specify positive and negative labels,
which is useful when you are sure the identity of a cell types depends
on the lack of expression of certain markers (see hematopoietic
lineages). To do so, you can add "+" or "-" to each gene.

# Downsampling

There are two ways to downsample data in `cellula`: downsampling
**reads** and downsampling **cells**.

The first approach simulates reads randomly sampling counts from a
distribution with a fixed total number, using a vector of probabilities
equivalent to the per-gene proportion of reads within each cell.

Briefly, let's consider a cell **C** in which genes *a*, *b*, and *c*
have been quantified with 50, 30, 20 counts each (totaling to 100
counts). This is equivalent to a bag of marbles in which the probability
of randomly picking an *a* marble is 50/100 = 0.5, *b* marble is 0.3,
and *c* marble is 0.2.

If we want to downsample **C** to a total of 40 counts (yielding the
downsampled **C'**), we randomly pick 40 counts from a (0.5, 0.3, 0.2)
vector of probabilities.

This is a sort of downsampling by simulation and is described in Scott
Tyler's [work](https://github.com/scottyler89/downsample)[[17]](#17),
reimplemented in `cellula` with a slightly faster optimization.

The `downsampleCounts()` uses a minimum count number that is
user-defined (or the minimum total count number in the dataset as a
default) and returns a `SingleCellExperiment` object with the same
number of cells as the input, and a down-sampled count matrix where each
cell has the same total number of counts.

The second approach randomly select cells from within groups such as
clusters, batches, or a combination of the two.

Cells are randomly selected so that they represent a user-defined
fraction of the within-group total, with some lower bound to ensure that
small groups are represented: if a rare cluster label only contains 9
cells and we want to downsample a dataset to 10%, we can cap the minimum
to 5 cells so that we ensure the rare label is still adequately
represented.

The `downsampleCells()` function returns a `SingleCellExperiment` object
with fewer cells than the input, as defined by the `proportion` and
`min` parameters.

# Inferring trajectories

At the time of writing `cellula` only implements a wrapper around the
`slingshot` method[[18]](#18) and the `monocle3` method [[19]](#19) for
pseudotemporal trajectory inference, and the `testPseudotime()` method
from `TSCAN`[[20]](#20) for differential expression along a trajectory.
More trajectory inference and differential expression methods will be
implemented in the future.

The `monocle3` method has been originally developed to work on 2D
embeddings such as UMAP, but given the distortion introduced by these
embeddings this function allows users to input any embedding (such as
PCA). Moreover, the `monocle3` method originally allows partitions to be
specified. In the current implementation partitions have been disabled.

The `findTrajectories()` function takes a `SingleCellExperiment` object
as input, and requires the user to specify the method (one of
`"slingshot"` or `"monocle"`), the cluster label which will be used as
an input to the MST creation in `slingshot`, or to identify the starting
node in `monocle3`.

Other important parameters are: 

  - `dr` - the reduced dimensional
reduction in which trajectories will be estimated (default is "PCA") 

  - `start` - the starting cluster for trajectory estimation, defaulting to
"auto" for an entropy-based method 

  - `omega` - `slingshot` only. Whether
or not to use a synthetic cluster to estimate disjointed trajectories,
as detailed in `?slingshot::getLineages`

  - `Monocle_lg_control` - `monocle` only. A list of control parameters 
  for the `learn_graph()` function from `monocle3`. 
  
  - `invert` - `monocle` only.
Whether or not to invert the direction of the pseudotime vector. 

  - `dr_embed` - an additional 2D dimensional reduction slot in which to
embed the trajectories (`slingshot` curves or `monocle` principal
graph). Yields slightly different results according to the method. 

  - `doDE` - whether or not to perform differential expression on every
trajectory.

The output is the same object used in the input, with some additional
fields.

-   the `colData` slot will contain a column containing the
    `slingPseudotime_N` pseudotime values for each lineage, where N is
    the number of the lineage (`slingshot`) or just a column named
    `monoclePseudotime` with pseudotime values for each cell
    (`monocle`).

-   the `metadata` slot will several new elements depending on the
    method.

    For `slingshot`:

    -   `Slingshot_lineages`: the output of `slingLineages()`, i.e. the
        result of the MST construction
    -   `Slingshot_embedded_curves`: the output of `embedCurves()`, i.e.
        the projection of low-dimensional principal curves onto a 2D
        embedding such as UMAP. Optional and only created if the
        parameter `dr_embed` is an available dimensionality reduction
        slot in the object.
    -   `Slingshot_MST`: the MST used to determine lineages by
        `slingshot`. Optional and only added if the parameter
        `add_metadata` is set to TRUE.
    -   `Slingshot_curves`: the principal curves constructed by
        `slingshot`. Optional and only added if the parameter
        `add_metadata` is set to TRUE.
    -   `Slingshot_weights`: the weights assigned per cell within each
        lineage by `slingshot`. Optional and only added if the parameter
        `add_metadata` is set to TRUE.
    -   `Slingshot_params`: the params used in the `slingshot` call.
        Optional and only added if the parameter `add_metadata` is set
        to TRUE.
    -   `pseudotime_DE`: the results of the differential expression
        test. This is a list where each lineage is tested separately.
        Optional and only created if `doDE` is set to TRUE.

    For `monocle`:

    -   `Monocle_embedded_curves`: the projection of the nodes of the
        principal graph to their nearest neighbor in 2D. Optional and
        only created if the parameter `dr_embed` is an available
        dimensionality reduction slot in the object.
    -   `Monocle_principal_graph`: the `igraph` object with the
        principal graph.
    -   `Monocle_principal_graph_aux`: the rest of the `CellDataSet`
        object with additional information on the principal graph.

In this example we download a dataset that captures early haematopoietic
differentiation from Nestorowa et al. [[21]](#21) and we do a quick
processing and clustering.

```{r}
sce2 = scRNAseq::NestorowaHSCData()

sce2 = cellula(sce2, name = "nestorowa", do_qc = FALSE)

sce2 = makeGraphsAndClusters(sce2, space = reducedDim(sce2, "PCA")[,1:20])

plot_UMAP(sce2, umap_slot = "UMAP", color_by = "SNN_0.64", label_by = "SNN_0.64")
```

<img src="https://github-production-user-asset-6210df.s3.amazonaws.com/21171362/242682308-dc8b13c4-f4ea-4dc1-93ff-e169ae1b64de.png" width="400"/>

Then we take an arbitrary starting point (cluster 7) and calculate
trajectories with both `slingshot` and `monocle` using the same space
(PCA, first 20 components).

```{r}
sce2 = findTrajectories(sce2, dr = "PCA", 
                        method = "monocle",
                        ndims = 20, clusters = "SNN_0.64", 
                        dr_embed = "UMAP", start = 7)
                        
sce2 = findTrajectories(sce2, dr = "PCA", 
                        method = "slingshot",
                        ndims = 20, clusters = "SNN_0.64", 
                        dr_embed = "UMAP", start = 7)           
```

Now we can visualize results on the UMAP since we specified a `dr_embed`
parameter. Notice how `slingshot` will identify different pseudotimes
(one per lineage) while `monocle3` will identify a single pseudotime.
Moreover, the `plot_UMAP()` function can include a `trajectory`
argument which points to the name of the element in the `metadata` slot
containing segments to draw trajectories.

```{r}
 plot_UMAP(sce2, umap_slot = "UMAP", color_by = "slingPseudotime_1", 
           label_by = "SNN_0.64", trajectory = "Slingshot_embedded_curves")
 plot_UMAP(sce2, umap_slot = "UMAP", color_by = "monoclePseudotime", 
           label_by = "SNN_0.64", trajectory = "Monocle_embedded_curves")
                                 
```

<img src="https://github-production-user-asset-6210df.s3.amazonaws.com/21171362/242682268-8855a6ca-c106-4f23-b2e6-c54fdcedac64.png" width=400>

<img src="https://github-production-user-asset-6210df.s3.amazonaws.com/21171362/242682289-8b43ed43-ab82-4d03-901a-d7109d3a3d52.png" width=400>

Since the 2D embedding of `monocle3` PCA-derived trajectories may be hard to 
understand, given the distortions introduced by UMAP, `cellula` includes an additional 
2D embedding method, `dr_embed = "FR"`, inspired båy the PAGA embedding 
initialization technique [[22]](#22).

Briefly, once the principal graph has been calculated, it is laid out in 2D 
using the Fruchterman-Reingold algorithm. Then each cell is randomly plotted 
around their closest vertex in the graph, and reordered according to pseudotime
value. This semi-random layout is used as an initialization for UMAP, which will
optimize the point positions. The resulting layout is more visually pleasing and
reflects more accurately the positions of cells with respect to the trajectory 
(although not necessarily to each other). 

This FR-initialized UMAP is stored in a `reducedDim` slot named `UMAP_FR`. 

```{r}
sce2 = findTrajectories(sce2, clusters = "SNN_0.64", method = "monocle",
        dr = "PCA", ndims = 20, start = "7", dr_embed = "FR")
        
plot_UMAP(sce2, umap_slot = "UMAP_FR", color_by = "monoclePseudotime", 
           label_by = "SNN_0.64", trajectories = "Monocle_embedded_curves")
           
```           
<img src="https://github-production-user-asset-6210df.s3.amazonaws.com/21171362/243069896-c1c933d0-8f23-4d73-b25c-5374eed73cb2.png" width=400>


It should be noted that this layout is, in a way, optimized for trajectories rather than global cell-cell similarity and as always should only be treated as a visualization tool.

## Metacells

In order to speed up calculations and overcome sparsity, cells can be
aggregated into metacells using k-means clustering with a high *k*.

Clustering is carried out in the reduced dimensional space of choice
selected through the `space` argument, and it is carried out in each
level of `group` separately.

Rather than selecting the number of clusters *k*, the function takes an
average number of cells per cluster *w* which is used to determine *k*
(default is w = 10 cells per cluster). Read counts are aggregated by
gene across all cells within each cluster, resulting in metacells. These
can be used as input to `findTrajectories()` or other operations such as
clustering, downsampling, signature scoring, etc.

In this example we create metacells aggregating (on average) 10 cells,
within each cluster from the "SNN_0.5" clustering results:

```{r}
sce_meta <- makeMetacells(sce, group = "SNN_0.5", space = "PCA_Harmony", w = 10)
```

# References

<a id="1">[1]</a> Amezquita et al. Nat Methods. 2020 Feb; 17(2):
137--145

<a id="2">[2]</a> Segerstolpe et al. Cell Metab. 2016 Oct 11; 24(4):
593--607.

<a id="3">[3]</a> Haghverdi et al. Nat Biotechnol. 2018 Jun;36(5):
421-427.

<a id="4">[4]</a> Lun ATL.
<https://marionilab.github.io/FurtherMNN2018/theory/description.html>

<a id="5">[5]</a> Korsunsky et al. Nat Methods. 2019 Dec;16(12):
1289-1296.

<a id="6">[6]</a> Stuart et al. Cell. 2019 Jun 13;177(7): 1888-1902.e21.

<a id="7">[7]</a> Welch et al. Cell. 2019 Jun 13; 177(7):
1873--1887.e17.

<a id="8">[8]</a> Germain et al. F1000Res. 2021; 10: 979

<a id="9">[9]</a> Blondel et al. J. Stat. Mech. 2008; P10008

<a id="10">[10]</a> Traag et al. Sci Rep. 2019; 9:5233.

<a id="11">[11]</a> Zappia and Oshlack Gigascience. 2018 Jul; 7(7):
giy083.

<a id="12">[12]</a> Aibar et al. Nat Methods. 2017 Nov; 14(11):
1083--1086.

<a id="13">[13]</a> Hänzelmann et al. BMC Bioinformatics. 2013 Jan
16;14: 7.

<a id="14">[14]</a> Andreatta and Carmona Comput Struct Biotechnol J.
2021 Jun 30;19: 3796-3798.

<a id="15">[15]</a> Jaitin et al. 2015 Science 343, 776-779

<a id="16">[16]</a> Muraro et al. Cell Syst. 2016 Oct 26;3(4):385-394.e3

<a id="17">[17]</a> Tyler et al. biorXiv 2021 11.15.468733

<a id="18">[18]</a> Street et al. BMC Genomics. 2018 Jun 19;19(1): 477.

<a id="19">[19]</a> Cao et al. Nature. 2019 Feb;566(7745):496-502

<a id="20">[20]</a> Ji and Ji
<https://www.bioconductor.org/packages/release/bioc/html/TSCAN.html>

<a id="21">[21]</a> Nestorowa et al. Blood. 2016 Aug 25;128(8):e20-31

<a id="22">[22]</a> Wolf et al. Genome Biol. 2019 Mar 19;20(1):59


