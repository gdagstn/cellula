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

    old_rownames = rownames(sce)
    rownames(sce) = paste0("g", seq_len(nrow(sce)))
    swap = FALSE

    if(any(duplicated(colnames(sce)))) {
      old_colnames = colnames(sce)
      swap = TRUE
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

    if(swap) colnames(sce) = old_colnames
    rownames(sce) = old_rownames

    if(verbose) cat(blue("[INT/Seurat]"), "Running UMAP on Seurat-corrected space.\n")

    reducedDim(sce, "UMAP_Seurat") <- umap(reducedDim(sce, "PCA_Seurat")[,seq_len(ndims)],
                                           n_neighbors = neighbor_n,
                                           min_dist = 0.7)
  } else if(method == "LIGER") {

    old_colnames = colnames(sce)

    if(any(duplicated(colnames(sce)))) {
      colnames(sce) = paste0("cell_", seq_len(ncol(sce)))
    }

    if(verbose) cat(blue("[INT/LIGER]"), "Creating LIGER object.\n")

    batches = unique(colData(sce)[,batch])
    countlist = lapply(batches, function(x) counts(sce[,colData(sce)[,batch] == x]))
    names(countlist) = batches
    l <- createLiger(countlist)
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
    reducedDim(sce, type = "LIGER") = do.call(rbind, l@H)
    reducedDim(sce, type = "LIGER_NORM") = l@H.norm

    if(verbose) cat(blue("[INT/LIGER]"), "Running UMAP on LIGER factorization.\n")
    reducedDim(sce, "UMAP_LIGER") <- umap(reducedDim(sce, "LIGER_NORM")[,seq_len(ndims)],
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
