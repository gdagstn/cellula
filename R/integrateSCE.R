#' SCE Integration pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels
#' @param hvgs character, vector of highly variable gene IDs
#' @param hvg_ntop numeric, the number of top HVGs to retain. Default is 2000.
#' @param method character, the integration method. One of \code{"fastMNN"}, \code{"Harmony"},
#'     \code{"Seurat"}, \code{"LIGER"}, and \code{"regression"}. Default is \code{"fastMNN"}.
#' @param ndims numeric, the number of dimensions to use for integration. Default
#'     is 20. For LIGER, it is the value of the k parameter in \code{rliger::optimizeALS()}
#' @param neighbor_n the number of neighbors used to compute UMAP. Default is
#'     NULL, which results in the rounded squared root of the number of cells.
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Note: for Seurat options, the
#'     \code{future} framework should be set up with maximum size and number of cores.
#' @param ret_nmf logical, should the NMF matrices be returned with the object? 
#'     if \code{TRUE}, they will be saved in the \code{metadata} of the object
#'     as a list named \code{"LIGER_NMF_matrices"}. Default is FALSE.
#' @param verbose logical, display messages on progress? Default is FALSE.
#'
#' @return a `SingleCellExperiment` object with integrated dimensionality reduction.
#'     In the case of LIGER integration, only the `H.norm` slot will be returned.
#'
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom scater runPCA
#' @importFrom uwot umap
#' @importFrom scran quickCluster computeSumFactors modelGeneVar
#' @importFrom batchelor multiBatchNorm fastMNN regressBatches
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom BiocParallel SerialParam
#'
#' @export

integrateSCE <- function(sce,
                        batch,
                        hvgs,
                        hvg_ntop = 2000,
                        method = "fastMNN",
                        ndims = 20,
                        neighbor_n = NULL,
                        parallel_param = SerialParam(),
                        ret_nmf = FALSE,
                        verbose = FALSE){
  ## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::integrateSCE()} - ")

  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if (!(method %in% c("fastMNN", "Harmony", "Seurat", "LIGER", "regression")))
    stop(ep, "method not recognized - must be one of \"fastMNN\", \"Harmony\", \"Seurat\", \"LIGER\", or \"regression\"")
  if (!(batch %in% colnames(colData(sce))))
    stop(ep,"batch column not found in the colData of the object")
  if (hvg_ntop > nrow(sce))
    stop(ep, "hvg_ntop cannot be higher than the number of features (nrow) in the object")

  if (is.null(neighbor_n)) {
    message(message(.bluem("[INT]"), "neighbor_n was not supplied. Defaulting to floor(sqrt(ncells))"))
    neighbor_n <- floor(sqrt(ncol(sce)))
  }

  # Start parameter logging - not fully implemented
  # TO DO
  # --------------------------------------------- #

  if (method == "fastMNN") {

    #fastMNN
    if (verbose) message(.bluem("[INT/fastMNN]"), "Correcting batch effect using fastMNN")
    if (is.null(hvgs)) {
      if (!is.null(metadata(sce)$hvgs)) {
        hvgs <- metadata(sce)$hvgs
        } else stop(ep, "Highly variable genes not supplied and not previously calculated")
      }

    sce_corr <- fastMNN(sce,
                        batch = colData(sce)[,batch],
                        d = ndims,
                        subset.row = hvgs,
                        BPPARAM = parallel_param)

    reducedDim(sce, "PCA_MNN") <- reducedDim(sce_corr, "corrected")


    # UMAP
    if (verbose) message(.bluem("[INT/fastMNN]"), "Running UMAP on MNN-corrected space")

    reducedDim(sce, "UMAP_MNN") <- umap(reducedDim(sce, "PCA_MNN")[,seq_len(ndims)],
                                        n_neighbors = neighbor_n,
                                        min_dist = 0.7)
  } else if (method == "Harmony") {


    if (!"harmony" %in% rownames(installed.packages()))
      stop(paste0(ep, "the `harmony` package must be installed first.\n
                Run `BiocManager::install(\"harmony\") to use this function."))

    if (verbose) message(.bluem("[INT/Harmony]"), "Correcting batch effect using Harmony")

    harmony_corr <- harmony::HarmonyMatrix(reducedDim(sce, "PCA")[,seq_len(ndims)],
                                                meta_data = colData(sce)[,batch],
                                                do_pca = FALSE)

    reducedDim(sce, "PCA_Harmony") <- harmony_corr

    if (verbose) message(.bluem("[INT/Harmony]"), "Running UMAP on Harmony-corrected space")

    reducedDim(sce, "UMAP_Harmony") <- umap(reducedDim(sce, "PCA_Harmony")[,seq_len(ndims)],
                                            n_neighbors = neighbor_n,
                                            min_dist = 0.7)


  } else if (method == "Seurat"){

    if (!"Seurat" %in% rownames(installed.packages()))
      stop(paste0(ep, "the `Seurat` package must be installed first.\n
                Run `BiocManager::install(\"Seurat\") to use this function."))

    if (verbose) message(.bluem("[INT/Seurat]"), "Converting to Seurat object.")
    old_colnames <- colnames(sce)

    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))

    nf <- data.frame(old_colnames, row.names = colnames(sce))

    seu <- SeuratObject::CreateSeuratObject(counts = counts(sce),
                                            meta.data = as.data.frame(colData(sce)))
    lc <- logcounts(sce)
    rownames(lc) <- rownames(seu)
    seu <- SeuratObject::SetAssayData(object = seu, slot = "data", new.data = lc)
    seu[["pca"]] <-SeuratObject:: CreateDimReducObject(embeddings = reducedDim(sce, "PCA"), key = "PC_")

    batches <-  unique(as.character((seu[[batch]][[batch]])))
    seulist <- lapply(batches, function(x) seu[,seu[[batch]] == x])
    names(seulist) <- batches

    if (verbose) message(.bluem("[INT/Seurat]"), "Normalization and HVG selection.")

    seulist <- lapply(seulist, function(x) {
      x <- Seurat::NormalizeData(x, verbose = verbose)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                        nfeatures = hvg_ntop,
                                        verbose = verbose)
    })

    if (verbose) message(.bluem("[INT/Seurat]"), "Finding anchors.")

    anchors <- Seurat::FindIntegrationAnchors(object.list = seulist,
                                              dims = seq_len(ndims),
                                              verbose = verbose,
                                              anchor.features = hvg_ntop)

    if (verbose) message(.bluem("[INT/Seurat]"), "Integration.")

    kweight <- min(100, unlist(lapply(seulist, function(x) floor(0.5*(ncol(x))))))

    if (verbose) message(.bluem("[INT/Seurat]"), " Using k.weight = ", kweight, ".")

    seu_int <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(ndims),
                             k.weight = kweight)

    if (verbose) message(.bluem("[INT/Seurat]"), "Running PCA on integrated object.")
    seu_int <- Seurat::ScaleData(seu_int,
                         assay = "integrated",
                         verbose = verbose)

    seu_int <- Seurat::RunPCA(seu_int,
                              npcs = ndims,
                              assay = "integrated",
                              verbose = verbose,
                              reduction.name = "spca")

    seu_int <- seu_int[,rownames(nf)]

    if (verbose) message(.bluem("[INT/Seurat]"), "Transferring to SCE object.")


    reducedDim(sce, "PCA_Seurat") <- SeuratObject::Embeddings(seu_int, reduction = "spca")[colnames(sce),]
    rm(seu)
    rm(seu_int)

    if (verbose) message(.bluem("[INT/Seurat]"), "Running UMAP on Seurat-corrected space.")

    reducedDim(sce, "UMAP_Seurat") <- umap(reducedDim(sce, "PCA_Seurat")[,seq_len(ndims)],
                                           n_neighbors = neighbor_n,
                                           min_dist = 0.7)
    colnames(sce) <- nf[colnames(sce), 1]

  } else if (method == "LIGER") {

    if (!"rliger" %in% rownames(installed.packages()))
      stop(paste0(ep, "the `rliger` and `RcppPlanc` packages must be installed first.\n
                Run `BiocManager::install(\"rliger\") to use this function."))

    old_colnames <- colnames(sce)

    if (any(duplicated(colnames(sce)))) {
      colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
    }

    nf <- data.frame(old_colnames, row.names = colnames(sce))

    if (verbose) message(.bluem("[INT/LIGER]"), "Creating LIGER object.")

    batches <- unique(colData(sce)[,batch])
    countlist <- lapply(batches, function(x) counts(sce[,colData(sce)[,batch] == x]))
    names(countlist) <- batches
    l <- rliger::createLiger(countlist)
    rm(countlist)

    if (verbose) message(.bluem("[INT/LIGER]"), "Preprocessing.")
    l <- rliger::normalize(l)
    l <- rliger::selectGenes(l)
    l <- rliger::scaleNotCenter(l)

    if (verbose) message(.bluem("[INT/LIGER]"), "Factorization.")
    infmats <- RcppPlanc::inmf(l@scale.data, k = ndims, verbose = verbose)
    l@H <- infmats$H
    l@W <- infmats$W
    l@V <- infmats$V
    
    for(i in seq_along(countlist)){
      rownames(l@H[[i]]) = colnames(l@scale.data)
      rownames(l@V[[i]]) = rownames(l@scale.data)
      rownames(l@W) = rownames(l@scale.data)
    }

    if (verbose) message(.bluem("[INT/LIGER]"), "Quantile normalization")
    l <- rliger::quantile_norm(l, verbose = verbose)

    if (verbose) message(.bluem("[INT/LIGER]"), "Transferring to SCE object.")
    reducedDim(sce, type = "LIGER") = do.call(rbind, l@H)[colnames(sce),]
    reducedDim(sce, type = "LIGER_NORM") = l@H.norm[colnames(sce),]
    
    if(ret_nmf) sce@metadata$LIGER_NMF_matrices = list("H" = l@H, "V" = l@V, "W" = l@W)

    if (verbose) message(.bluem("[INT/LIGER]"), "Running UMAP on LIGER factorization.")
    reducedDim(sce, "UMAP_LIGER") <- umap(reducedDim(sce, "LIGER_NORM")[,seq_len(min(ncol(reducedDim(sce, "LIGER_NORM")), ndims))],
                                          n_neighbors = neighbor_n,
                                          min_dist = 0.7)

    colnames(sce) <- nf[colnames(sce), 1]

  } else if (method == "regression") {

    #regression
    if (verbose) message(.bluem("[INT/regression]"), "Correcting batch effect using regression")

    sce_corr <- regressBatches(sce,
                               batch = colData(sce)[,batch],
                               subset.row = hvgs,
                               d = ndims,
                               BPPARAM = parallel_param)

    reducedDim(sce, "PCA_regression") <- reducedDim(sce_corr, "corrected")

    if (verbose) message(.bluem("[INT/regression]"), "Running UMAP on regression-corrected space.")

    reducedDim(sce, "UMAP_regression") <- umap(reducedDim(sce, "LIGER_NORM")[,seq_len(min(ncol(reducedDim(sce, "PCA_regression")), ndims))],
                                          n_neighbors = neighbor_n,
                                          min_dist = 0.7)
  }
  return(sce)
}
