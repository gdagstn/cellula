#'#' SCE Integration pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels
#' @param name character, the name of the project/folder where files will be saved.
#'     If NULL, a random name will be generated
#' @param geneset_list named list of gene IDs (must coincide with rownames of the
#'     SCE object) that will be used by AUCell. Default is NULL meaning no AUC
#'     will be calculated.
#' @param discard logical, should values that do not meet QC thresholds be discarded?
#'     Default is TRUE.
#' @param subset_mito logical, should mitochondrial transcripts be used for QC?
#'     Default is TRUE.
#' @param subset_ribo logical, should ribosomal transcripts be used for QC?
#'     Default is TRUE.
#' @param subset_malat1 logical, should MALAT1 transcripts be used for QC?
#'     Default is TRUE.
#' @param detect_doublets logical, should `scDblFinder` be run? Default is TRUE.
#' @param run_emptydrops logical, should `emptyDrops` be run? Default is TRUE.
#' @param emptydrops_cutoff either "auto" (default, barcode rank inflection point)
#'     or a numeric. Cells with total reads below this cutoff are used to calculate
#'     ambient RNA profiles and are removed.
#' @param emptydrops_alpha numeric, the FDR threshold to call an empty barcode.
#' @param hvg_ntop numeric, the number of top highly variable genes to be selected.
#'     Default is 2000.
#' @param integration_method character, one of "fastMNN", "Harmony", "Seurat",
#'     "LIGER", or "regression".
#' @param ndims numeric, the number of dimensions to retain in the reduced dimension
#'     embedding for downstream applications. Default is 20.
#' @param verbose logical, display messages on progress? Default is FALSE.
#' @param save_plots logical, should plots be drawn and saved? Default is TRUE
#' @param parallel_param a BiocParallel object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Note: for Seurat options, the
#'     `future` framework should be set up with maximum size and number of cores.
#'
#' @return  a `SingleCellExperiment` object with normalized data, doublet assignment
#'    (if calculated), uncorrected and corrected PCA and UMAP coordinates according
#'    to the method of choice.
#'
#' @importFrom ids adjective_animal
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom DropletUtils emptyDrops barcodeRanks
#' @importFrom crayon blue
#' @importFrom scuttle isOutlier perCellQCMetrics
#' @importFrom gridExtra grid.arrange
#' @importFrom scater plotColData runPCA plotReducedDim
#' @importFrom uwot umap
#' @importFrom ggplot2 scale_y_log10 ggtitle ggsave
#' @importFrom scDblFinder scDblFinder
#' @importFrom scran quickCluster computeSumFactors modelGeneVar getTopHVGs
#' @importFrom batchelor multiBatchNorm
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom BiocParallel SerialParam
#'
#' @export

pipeline_integration <- function(sce,
                                 batch,
                                 name = NULL,
                                 geneset_list = NULL,
                                 discard = TRUE,
                                 subset_mito = TRUE,
                                 subset_ribo = TRUE,
                                 subset_malat1 = TRUE,
                                 detect_doublets = TRUE,
                                 run_emptydrops = FALSE,
                                 emptydrops_cutoff = "auto",
                                 emptydrops_alpha = 0.01,
                                 hvg_ntop = 2000,
                                 integration_method = "fastMNN",
                                 ndims = 20,
                                 verbose = FALSE,
                                 save_plots = TRUE,
                                 parallel_param = SerialParam()) {

  if(!batch %in% colnames(colData(sce)))
    stop(paste0("Batch label \"", batch, "\" not found."))

  # Make folder
  if(is.null(name)) {
    name = adjective_animal()
    warning(paste0("No name selected so your folder name will be: ", name))
  }

  dir.create(name)

  if(!is.factor(colData(sce)[,batch]))
    colData(sce)[,batch] = as.factor(colData(sce)[,batch])

  if(verbose) {
    cat("Working on object", name, "\n")
    ncells = ncol(sce)
    cat("Input cells: ", ncells, "\n")
    cat("By batch:\n")
    print(table(colData(sce)[,batch]))
  }

  if(run_emptydrops){

    if(verbose) cat(blue("[EMPTY]"),"Running emptyDrops. \n")

    if(emptydrops_cutoff == "auto") {
      barcode_ranks <- barcodeRanks(sce)
      emptydrops_cutoff = metadata(barcode_ranks)$inflection
    }

    empty_droplets <- emptyDrops(sce, lower = emptydrops_cutoff)
    keep_droplets <- empty_droplets$FDR <= emptydrops_alpha
    sce$empty <- factor(ifelse(empty_droplets$FDR <= emptydrops_alpha, "ok", "empty"))
    sce$empty[which(is.na(sce$empty))] = "empty"
    if(verbose) print(table(Sig = keep_droplets, Limited = empty_droplets$Limited))
    if(verbose) cat("Empty cells: ", sum(sce$empty == "empty"), "\n")
    sce = sce[, which(sce$empty == "ok"), drop = FALSE]
  }

  if(verbose) cat(blue("[QC]"),"Calculating QC metrics. \n")

  # Generate QC
  if(subset_mito){
    if(sum(grepl("^MT", rowData(sce)$Symbol)) == 0) {
      cat("   No MT genes found.\n")
      mito = NA
    } else mito = rownames(sce)[grepl("^MT-", rowData(sce)$Symbol)]
  } else mito = NA

  if(subset_malat1){
    if(sum(grepl("^MALAT1", rowData(sce)$Symbol)) == 0) {
      cat("   No MALAT1 gene found.\n")
      Malat1 = NA
    } else Malat1 = rownames(sce)[grepl("^MALAT1", rowData(sce)$Symbol)]
  } else Malat1 = NA

  if(subset_ribo){
    if(sum(grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol)) == 0) {
      cat("   No ribo genes found.\n")
      Ribo = NA
    } else Ribo = rownames(sce)[grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol)]
  } else Ribo = NA

  subset_list = list(mito = mito, Malat1 = Malat1, Ribo = Ribo)
  subset_list = subset_list[!is.na(subset_list)]

  if(length(subset_list) > 0) {
    sce_fqc <- perCellQCMetrics(sce,
                                subsets = subset_list,
                                BPPARAM = parallel_param)

    low.lib <- isOutlier(log10(sce_fqc$sum),
                         batch = colData(sce)[,batch],
                         type = "lower",
                         nmads=3)

    low.genes <- isOutlier(log10(sce_fqc$detected),
                           batch = colData(sce)[,batch],
                           type = "lower",
                           nmads=3)

    if(!all(is.na(mito))) {
      high.mt <- isOutlier(sce_fqc$subsets_mito_percent,
                           batch = colData(sce)[,batch],
                           type = "higher",
                           nmads = 3)
    } else high.mt <- NA

    if(!all(is.na(Malat1))) {
      high.malat1 <- isOutlier(sce_fqc$subsets_Malat1_percent,
                               batch = colData(sce)[,batch],
                               type = "higher",
                               nmads = 3)
    } else high.malat1 <- NA

    if(!all(is.na(Ribo))) {
      high.ribo <- isOutlier(sce_fqc$subsets_Ribo_percent,
                             batch = colData(sce)[,batch],
                             type = "higher",
                             nmads = 3)
    } else high.ribo <- NA

    data.frame(LowLib=sum(low.lib),
               LowNgenes=sum(low.genes),
               HighMT = ifelse(all(is.na(mito)), FALSE, sum(high.mt)),
               HighMalat1 = ifelse(all(is.na(Malat1)), FALSE, sum(high.malat1)),
               HighRibo = ifelse(all(is.na(Ribo)), FALSE, sum(high.ribo)))

    sce_fqc$discard <- low.lib | low.genes | high.mt | high.malat1 | high.ribo

    colData(sce) <- cbind(colData(sce), sce_fqc)

    if(save_plots){

      if(verbose) cat(blue("[QC]"),"   Saving QC plots. \n")

      p1 <- grid.arrange(
        plotColData(sce, y="sum", x = batch, colour_by="discard") +
          scale_y_log10() + ggtitle("Total count"),
        plotColData(sce, y="detected", x = batch, colour_by="discard") +
          scale_y_log10() + ggtitle("Detected features"),
        plotColData(sce, y="subsets_Malat1_percent", x = batch,
                    colour_by="discard") + ggtitle("Malat1 percent"),
        plotColData(sce, y="subsets_Ribo_percent", x = batch,
                    colour_by="discard") + ggtitle("Ribo percent"),
        plotColData(sce, y="subsets_mito_percent", x = batch,
                    colour_by="discard") + ggtitle("Mito percent"),
        ncol=2
      )
      ggsave(p1, filename = paste0("./", name, "/QC_plot.png"), width =  6, height = 6, device = "png")
    }

    sce$discard[is.na(sce$discard)] = TRUE

    if(discard) sce <- sce[, !sce$discard]
  }

  # Doublet finding
  if(detect_doublets){
    if(verbose) cat(blue("[DBL]"), "Finding doublets. \n")

    sce <- scDblFinder(sce, verbose = verbose,
                       samples = colData(sce)[,batch],
                       BPPARAM = parallel_param)
    if(save_plots) {
      p2 <- plotColData(sce, x="scDblFinder.class", y = "sum", colour_by="scDblFinder.class") +
        scale_y_log10() + ggtitle("Doublets")
      ggsave(p2, filename = paste0("./", name, "/doublet_plot.png"), width =  6, height = 6, device = "png")
    }
  }
  #if(!is.null(stopat) & stopat == "DBL") return(sce)
  if(verbose) cat("Saving temporary file. \n")

  saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))

  if(verbose) cat(blue("[NORM]"), "Calculating size factors and normalizing. \n")

  # Size factors
  if(verbose) cat(blue("[NORM]"),"   Preclustering. \n")

  sce_cl <- quickCluster(sce,
                         block = colData(sce)[,batch],
                         BPPARAM = parallel_param)

  if(verbose) cat(blue("[NORM]"),"   Calculating pooled factors. \n")
  sce <- computeSumFactors(sce,
                           clusters = sce_cl,
                           BPPARAM = parallel_param)

  # Normalization
  if(verbose) cat(blue("[NORM]"),"   Log-normalization. \n")

  sce <- multiBatchNorm(sce,
                        batch = colData(sce)[,batch])

  if(verbose) cat("Saving temporary file. \n")

  saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))

  #if(!is.null(stopat) & stopat == "NORM") return(sce)

  # HVGs

  if(verbose) cat(blue("[DR]"), "Selecting HVGs. \n")

  vargenes <- modelGeneVar(sce, block = colData(sce)[,batch])

  hvgs <- getTopHVGs(vargenes, n = hvg_ntop)

  # PCA

  if(verbose) cat(blue("[DR]"), "Running PCA. \n")

  sce <- runPCA(sce,
                subset_row = hvgs,
                exprs_values = "logcounts")#,
  #BPPARAM = parallel_param)

  if(verbose) cat("Saving temporary file. \n")

  saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))

  # UMAP
  if(verbose) cat(blue("[DR]"), "Running UMAP on uncorrected PCA \n")

  neighbor_n <- sqrt(ncol(sce))

  reducedDim(sce, "UMAP") <- umap(reducedDim(sce, "PCA")[,seq_len(ndims)],
                                                               n_neighbors = neighbor_n,
                                                               min_dist = 0.7)

  if(verbose) cat("Saving temporary file. \n")

  saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))

  if(verbose) cat(blue("[INT]"), "Integration. \n")

  sce = integrateSCE(sce,
                     batch = batch,
                     hvgs = hvgs,
                     hvg_ntop = hvg_ntop,
                     method = integration_method,
                     ndims = ndims,
                     parallel_param = parallel_param,
                     verbose = verbose)

  if(verbose) cat("Done.\nSaving temporary file. \n")

  saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))

  # AUCell
  if(!is.null(geneset_list)) {

    if(verbose) cat("Assigning cell labels. \n")

    rankings <- AUCell_buildRankings(counts(sce),
                                     plotStats = FALSE,
                                     verbose = FALSE)

    aucs <- AUCell_calcAUC(geneset_list,
                           rankings,
                           aucMaxRank = ceiling(0.2 * nrow(rankings)))
    # All assignments

    assigned <- as.data.frame(t(assay(aucs)))


    # Best overall score
    assigned$first_max_score <- apply(assigned, 1, max)

    # Second best score
    assigned$second_max_score <- apply(assigned[,1:(ncol(assigned) - 1)], 1, function(x) {
      return(max(x[x != max(x)]))
    })

    # Ambiguous labels
    assigned$ambiguous <- (assigned$first_max_score - assigned$second_max_score)/(assigned$first_max_score + assigned$second_max_score) <= 0.2

    assigned$best_label = colnames(assigned)[apply(assigned[,1:(ncol(assigned) - 3)], 1, which.max)]

    colData(sce)$labels <- factor(assigned$best_label)


    p3 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "labels") + ggtitle("Labels")

    ggsave(p3, filename = paste0("./", name, "/UMAP_labels_plot.png"), width =  6, height = 6, device = "png")

  }
  if(verbose) cat("Saving final object.\n")
  saveRDS(sce, file = paste0("./", name, "/", name, "_PS_INT_SCE.RDS"))

  if(verbose) cat("Deleting temporary file. \n")
  file.remove(paste0("./", name, "/", name, "_tempSCE.RDS"))

  if(verbose) cat("All done. Input cells: ", ncells, ", final cell number: ", ncol(sce), ".\n")
  return(sce)
}


#' SCE multi-resolution clustering
#'
#' Very simple wrapper to SNN graph and Louvain/Leiden clustering using
#' multiple resolutions
#'
#' @param sce a SingleCellExperiment object
#' @param neighbors numeric, number of neighbors for SNN graph edge construction.
#'     Default is 10.
#' @param weighting_scheme character, the weighting scheme for SNN graph construction.
#'     One of "jaccard", "rank", "number". Default is "jaccard".
#' @param sweep_on character, the parameter used for sweeping. Can be "clustering",
#'     meaning values of `k` will be looped through as resolution, or "SNN",
#'     meaning values of `k`will be looped through as number of neighbors.
#' @param method character, the type of graph-based clustering to use. One of
#'     "louvain" or leiden". Default is "louvain".
#' @param k numeric, vector of parameter sweep for graph cosntruction or clustering.
#' @param space a matrix of lower dimensional embedding such as the PCA coordinates.
#'     if NULL (default), the "PCA" slot from `reducedDims(sce)`
#' @param ndims numeric, the number of dimensions (columns of `space`) to use to
#'     build the SNN graph. Default is 20
#' @param calculate_modularity logical, should pairwise modularity between
#'     clusters be calculated? Default is TRUE
#' @param calculate_silhouette logical, should approximate silhouette widths be
#'     calculated? Default is TRUE
#' @param leiden_iterations numeric, the number of iterations of Leiden clustering.
#'     Default is 5.
#' @param prefix character, the prefix of the column names on `colData`
#'     where clustering results are stored. Default is "SNN_".
#' @param verbose logical, should messages be written? Default is FALSE
#' @param BPPARAM a `BiocParallel` parameter object for graph construction
#'
#' @return a `SingleCellExperiment` object with cluster memberships in the
#'     `colData` table, named according to prefix and respective value of `k`.
#'     Optionally, silhouette and/or modularity values are stored in the `metadata`
#'     slot of the `SingleCellExperiment` object, one for every value of `k`.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom crayon blue
#' @importFrom bluster makeSNNGraph pairwiseModularity approxSilhouette
#' @importFrom igraph cluster_louvain cluster_leiden
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom BiocParallel SerialParam
#'
#'
#' @export

makeGraphsAndClusters <- function(sce,
                                  neighbors = 10,
                                  weighting_scheme = "jaccard",
                                  sweep_on = "clustering",
                                  method = "louvain",
                                  k = seq(0.1, 1, length.out = 6),
                                  space = NULL,
                                  ndims = 20,
                                  calculate_modularity = TRUE,
                                  calculate_silhouette = TRUE,
                                  leiden_iterations = 5,
                                  prefix = "SNN_",
                                  verbose = FALSE,
                                  BPPARAM = SerialParam()) {

  if(is.null(space)) space = reducedDim(sce, "PCA")[,seq_len(ndims)]


  # Case 1: parameter sweep on clustering resolution

  if(sweep_on == "clustering" & !is.null(neighbors)) {

    if(verbose) cat(blue("[CLU]"), "Creating SNN graph.\n")

    g = makeSNNGraph(space,
                     k = neighbors,
                     type = weighting_scheme,
                     BPPARAM = BPPARAM)
  for(i in k){

     if(verbose) cat(blue("[CLU]"), "Clustering at resolution ", i, ".\n")

     if(method == "louvain") {
      cl = factor(cluster_louvain(g, resolution = i)$membership)
     } else if(method == "leiden") {
      cl = factor(cluster_leiden(g, objective_function = "CPM",
                                         n_iterations = leiden_iterations,
                                         resolution_parameter = i)$membership)
     }

    gname = paste0(prefix, i)
    colData(sce)[,gname] = cl

    if(verbose) cat(blue("[CLU]"), "Found", length(unique(cl)), "clusters.\n")

    if(calculate_modularity) {
      if(verbose) cat(blue("[CLU]"), "Calculating pairwise modularity.\n")
      metadata(sce)[[paste0("modularity_", gname)]] = pairwiseModularity(g, clusters = cl, as.ratio = TRUE)
    }

    if(calculate_silhouette) {
      if(verbose) cat(blue("[CLU]"), "Calculating approximate silhouette widths.\n")
      silhouette <- as.data.frame(approxSilhouette(space, clusters = cl))
      silhouette$closest <- factor(ifelse(silhouette$width > 0, cl, silhouette$other))
      silhouette$cluster <- cl
      metadata(sce)[[paste0("silhouette_", gname)]] = silhouette
        }
    }
    # Case 2: parameter sweep on SNN neighbor number
  } else if(sweep_on == "SNN" & !is.null(k)) {

    k = floor(k)

    for(i in k) {
      if(verbose) cat(blue("[CLU]"), "Creating SNN graph with k =", i, "neighbors.\n")
      g = makeSNNGraph(space,
                       k = i,
                       type = weighting_scheme,
                       BPPARAM = BPPARAM)
      if(method == "louvain") {
        cl = factor(cluster_louvain(g)$membership)
      } else if(method == "leiden") {
        cl = factor(cluster_leiden(g, objective_function = "CPM",
                                   n_iterations = leiden_iterations)$membership)
      }
      gname = paste0(prefix, i)
      colData(sce)[,gname] = cl
      if(verbose) cat("Found", length(unique(cl)), "clusters.\n")

      if(calculate_modularity) {
        if(verbose) cat(blue("[CLU]"), "Calculating pairwise modularity.\n")
        metadata(sce)[[paste0("modularity_", gname)]] = pairwiseModularity(g, clusters = cl, as.ratio = TRUE)
      }

      if(calculate_silhouette) {
        if(verbose) cat(blue("[CLU]"), "Calculating approximate silhouette widths.\n")
        sil <- as.data.frame(approxSilhouette(space, clusters = cl))
        sil$closest <- factor(ifelse(sil$width > 0, cl, sil$other))
        sil$cluster <- cl
        metadata(sce)[[paste0("silhouette_", gname)]] = sil
      }
     }
  }

 return(sce)
}



#' Plot cluster modularity
#'
#' Plots a heatmap showing pairwise cluster modularity from a SCE object
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "modularity_SNN_100"
#'
#' @return a heatmap of pairwise modularity
#'
#' @importFrom pheatmap pheatmap
#' @importFrom colorspace sequential_hcl
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export


plotModularity <- function(sce, name) {
  name = paste0("modularity_", name)
  pheatmap(mat = log2(metadata(sce)[[name]] + 1),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = sequential_hcl(palette = "Sunset", n = 25),
                     main = "log2(modularity ratio + 1)",
                     border_color = NA
                    )
}


#' Plot approximate silhouette widths
#'
#' Plots a beesswarm plot showing the approximate width of each cell within
#' each cluster
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "silhouette_SNN_100"
#'
#' @return a beeswarm plot of silhouette widhts
#'
#' @importFrom ggplot2 ggplot aes_string theme_bw
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export

plotSilhouette <- function(sce, name) {
 name = paste0("silhouette_", name)
 ggplot(metadata(sce)[[name]], aes_string(x="cluster", y="width", colour="closest")) +
          geom_quasirandom(method="smiley") +
    theme_bw()
}

#' SCE Integration pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels
#' @param hvgs character, vector of highly variable gene IDs
#' @param hvg_ntop numeric, the number of top HVGs to retain. Default is 2000.
#' @param method character, the integration method. One of "fastMNN", "Harmony",
#'     "Seurat", "LIGER", and "regression". Default is "fastMNN".
#' @param ndims numeric, the number of dimensions to use for integration. Default
#'     is 20.
#' @param liger_k the k parameter for LIGER integration, i.e. dimensionality of
#'     the cell resulting integrated embeddings. Default is 20.
#' @param neighbor_n the number of neighbors used to compute UMAP. Default is
#'     NULL, which results in the rounded squared root of the number of cells.
#' @param verbose logical, display messages on progress? Default is FALSE.
#' @param parallel_param a BiocParallel object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Note: for Seurat options, the
#'     `future` framework should be set up with maximum size and number of cores.
#'
#' @return a `SingleCellExperiment` object with integrated dimensionality reduction.
#'     In the case of LIGER integration, only the `H.norm` slot will be returned.
#'
#' @importFrom ids adjective_animal
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom crayon blue
#' @importFrom scater plotColData runPCA
#' @importFrom uwot umap
#' @importFrom scDblFinder scDblFinder
#' @importFrom scran quickCluster computeSumFactors modelGeneVar
#' @importFrom batchelor multiBatchNorm fastMNN regressBatches
#' @importFrom harmony HarmonyMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData NormalizeData
#'     FindVariableFeatures FindIntegrationAnchors IntegrateData ScaleData RunPCA
#' @importFrom SeuratObject CreateDimReducObject Embeddings
#' @importFrom rliger createLiger normalize selectGenes scaleNotCenter optimizeALS
#'     quantile_norm
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom BiocParallel SerialParam
#'
#' @export

integrateSCE = function(sce,
                        batch,
                        hvgs,
                        hvg_ntop = 2000,
                        method = "fastMNN",
                        ndims = 20,
                        liger_k = 20,
                        neighbor_n = NULL,
                        parallel_param = SerialParam(),
                        verbose = FALSE){

  if(is.null(neighbor_n)) neighbor_n = floor(sqrt(ncol(sce)))

  if(method == "fastMNN") {

    #fastMNN
    if(verbose) cat(blue("[INT/fastMNN]"), "Correcting batch effect using fastMNN.\n")

    sce_corr <- fastMNN(sce,
                        batch = colData(sce)[,batch],
                        subset.row = hvgs,
                        BPPARAM = parallel_param)

    reducedDim(sce, "PCA_MNN") <- reducedDim(sce_corr, "corrected")


    # UMAP
    if(verbose) cat(blue("[INT/fastMNN]"), "Running UMAP on MNN-corrected space.\n")

    reducedDim(sce, "UMAP_MNN") <- umap(reducedDim(sce, "PCA_MNN")[,seq_len(ndims)],
                                        n_neighbors = neighbor_n,
                                        min_dist = 0.7)
  } else if(method == "Harmony") {

    if(verbose) cat(blue("[INT/Harmony]"), "Correcting batch effect using Harmony.\n")

    harmony_corr <- HarmonyMatrix(reducedDim(sce, "PCA")[,seq_len(ndims)],
                                  meta_data = colData(sce)[,batch],
                                  do_pca = FALSE)

    reducedDim(sce, "PCA_Harmony") <- harmony_corr

    if(verbose) cat(blue("[INT/Harmony]"), "Running UMAP on Harmony-corrected space.\n")

    reducedDim(sce, "UMAP_Harmony") <- umap(reducedDim(sce, "PCA_Harmony")[,seq_len(ndims)],
                                            n_neighbors = neighbor_n,
                                            min_dist = 0.7)


  } else if(method == "Seurat"){

    if(verbose) cat(blue("[INT/Seurat]"), "Converting to Seurat object.\n")

    if(any(duplicated(colnames(sce)))) {
      old_colnames = colnames(sce)
      colnames(sce) = paste0("cell_", seq_len(ncol(sce)))
    }
    seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
    seu <- SetAssayData(object = seu, slot = "data", new.data = logcounts(sce))
    seu[["pca"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "PCA"), key = "PC_")

    batches =  unique(as.character((seu[[batch]][[batch]])))
    seulist = lapply(batches, function(x) seu[,seu[[batch]] == x])
    names(seulist) = batches

    if(verbose) cat(blue("[INT/Seurat]"), "Normalization and HVG selection.\n")

    seulist = lapply(seulist, function(x) {
      x <- NormalizeData(x, verbose = verbose)
      x <- FindVariableFeatures(x,
                                selection.method = "vst",
                                nfeatures = hvg_ntop,
                                verbose = verbose)
    })

    if(verbose) cat(blue("[INT/Seurat]"), "Finding anchors.\n")

    anchors <- FindIntegrationAnchors(object.list = seulist,
                                      dims = seq_len(ndims),
                                      verbose = verbose)

    if(verbose) cat(blue("[INT/Seurat]"), "Integration.\n")

    seu_int <- IntegrateData(anchorset = anchors, dims = seq_len(ndims))

    if(verbose) cat(blue("[INT/Seurat]"), "Running PCA on integrated object.\n")
    seu_int <- ScaleData(seu_int,
                         assay = "integrated",
                         verbose = verbose)

    seu_int <- RunPCA(seu_int,
                      npcs = ndims,
                      assay = "integrated",
                      verbose = verbose,
                      reduction.name = "spca")

    if(verbose) cat(blue("[INT/Seurat]"), "Transferring to SCE object.\n")

    reducedDim(sce, "PCA_Seurat") = Embeddings(seu_int, reduction = "spca")
    rm(seu)
    rm(seu_int)

    if(verbose) cat(blue("[INT/Seurat]"), "Running UMAP on Seurat-corrected space.\n")

  reducedDim(sce, "UMAP_Seurat") <- umap(reducedDim(sce, "PCA_Seurat")[,seq_len(ndims)],
                                         n_neighbors = neighbor_n,
                                         min_dist = 0.7)
    colnames(sce) = old_colnames
  } else if(method == "LIGER") {

    if(any(duplicated(colnames(sce)))) {
      old_colnames = colnames(sce)
      colnames(sce) = paste0("cell_", seq_len(ncol(sce)))
    }

    if(verbose) cat(blue("[INT/LIGER]"), "Creating LIGER object.\n")

    batches = unique(colData(sce)[,batch])
    countlist = lapply(batches, function(x) counts(sce[,colData(sce)[,batch] == x]))
    names(countlist) = batches
    l <- (countlist)
    rm(countlist)

    if(verbose) cat(blue("[INT/LIGER]"), "Preprocessing.\n")
    l <- normalize(l)
    l <- selectGenes(l)
    l <- scaleNotCenter(l)

    if(verbose) cat(blue("[INT/LIGER]"), "Factorization.\n")
    l <- optimizeALS(l, k = liger_k, verbose = verbose)

    if(verbose) cat(blue("[INT/LIGER]"), "Quantile normalization\n")
    l <- quantile_norm(l, verbose = verbose)

    if(verbose) cat(blue("[INT/LIGER]"), "Transferring to SCE object.\n")
    reducedDim(sce, type = "LIGER") = l@H.norm

    if(verbose) cat(blue("[LIGER]"), "Running UMAP on LIGER factorization.\n")
    reducedDim(sce, "UMAP_LIGER") <- umap(reducedDim(sce, "LIGER")[,seq_len(ndims)],
                                          n_neighbors = neighbor_n,
                                          min_dist = 0.7)

    colnames(sce) = old_colnames

  } else if(method == "regression") {

    #regression
    if(verbose) cat(blue("[INT/regression]"), "Correcting batch effect using regression\n")

    sce_corr <- regressBatches(sce,
                               batch = colData(sce)[,batch],
                               subset.row = hvgs,
                               d = ndims,
                               BPPARAM = parallel_param)

    reducedDim(sce, "PCA_regression") <- reducedDim(sce_corr, "corrected")

  }
  return(sce)
}
