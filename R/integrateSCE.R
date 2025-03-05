#' SCE Integration pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels
#' @param hvgs character, vector of highly variable gene IDs
#' @param hvg_ntop numeric, the number of top HVGs to retain. Default is 2000.
#' @param method character, the integration method. One of \code{"fastMNN"}, \code{"Harmony"},
#'     \code{"Seurat"}, \code{"LIGER"}, \code{"regression"}, \code{"scMerge2"},
#'     and \code{"STACAS"}. Default is \code{"fastMNN"}.
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
#' @param ... extra arguments passed to the main integration functions used by each method.
#'
#' @return a `SingleCellExperiment` object with integrated dimensionality reduction.
#'     In the case of LIGER integration, only the "H" and "H.norm" slots will be returned,
#'     unless the option \code{ret_nmf} is set to \code{TRUE}, in which case the
#'     `H`, `W` and `V` matrices will be returned as well in the \code{metadata} slot.
#'     The integrated dimensionality reductions will be called "UMAP_ or PCA_\{method\}"
#'     where \{method\} is any method other than "LIGER", in which case the
#'     dimensionality reduction are called "LIGER", "LIGER_NORM" and "UMAP_LIGER".
#'
#' @details
#' This function allows to use several dataset integration/batch effect
#' correction methods for single cell datasets, allowing the user to supply
#' any number of parameters using \code{...}. Due to syntax limitations, the
#' parameters that can be tuned will only apply to *one* function within any
#' integration pipeline. In particular, the following functions are affected:
#'
#'\itemize{
#'  \item{method = "fastMNN":}{ \code{batchelor::fastMNN()}}
#'  \item{method = "Harmony":}{ \code{harmony::RunHarmony()}}
#'  \item{method = "Seurat":}{  \code{Seurat::FindIntegrationAnchors()}}
#'  \item{method = "LIGER":}{ \code{RcppPlanc::inmf}}
#'  \item{method = "regression":}{ \code{batchelor::regressBatches()}}
#'  \item{method = "scMerge2"}{ \code{scMerge::scMerge2()}}
#'  \item{method = "STACAS"}{ \code{STACAS::RunStacas()}}
#'}
#'
#' Users who desire further control/customization should apply the functions
#' from the respective packages directly, e.g. to change several parameters
#' in the \code{Seurat} integration pipeline.
#'
#' The following methods work out of the box:
#'
#'\itemize{
#'  \item{method = "fastMNN":}{ FastMNN correction from \code{batchelor}}
#'  \item{method = "Harmony":}{Integration on PCA embeddings from \code{harmony}}
#'  \item{method = "Seurat":}{\code{Seurat} CCA with de novo normalization
#'      and feature selection, anchor finding and integration}
#'  \item{method = "LIGER":}{LIGER iNMF using the \code{RcppPlanc} implementation,
#'      with de novo normalization and feature selection through \code{rliger}}
#'  \item{method = "regression":}{Linear regression from \code{batchelor}}
#'  \item{method = "scMerge2"}{scMerge2 pseudobulking and RUV from \code{scMerge}}
#'  \item{method = "STACAS"}{\code{Seurat} pre-processing and \code{STACAS} integration}
#'}
#'
#' @importFrom SummarizedExperiment colData
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
                        verbose = FALSE,
                        ...){
  ## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::integrateSCE()} - ")

  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if (!(method %in% c("fastMNN", "Harmony", "Seurat", "LIGER", "regression", "STACAS", "scMerge2")))
    stop(ep, "method not recognized - must be one of \"fastMNN\", \"Harmony\",
         \"Seurat\", \"LIGER\", \"regression\", \"scMerge2\", \"STACAS\"")
  if (!(batch %in% colnames(colData(sce))))
    stop(ep,"batch column not found in the colData of the object")
  if (hvg_ntop > nrow(sce))
    stop(ep, "hvg_ntop cannot be higher than the number of features (nrow) in the object")

  if (is.null(neighbor_n)) {
    message(message(.bluem("[INT]"), " neighbor_n was not supplied. Defaulting to floor(sqrt(ncells))"))
    neighbor_n <- floor(sqrt(ncol(sce)))
  }

  # Start parameter logging - not fully implemented
  # TO DO
  # --------------------------------------------- #

  if (method == "fastMNN") {

    sce <- .integrateMNN(sce = sce,
                         batch = batch,
                         hvgs = hvgs,
                         ndims = ndims,
                         neighbor_n = neighbor_n,
                         parallel_param = parallel_param,
                         verbose = verbose,
                         ...)

  } else if (method == "Harmony") {

    sce <- .integrateHarmony(sce = sce,
                             batch = batch,
                             ndims = ndims,
                             neighbor_n = neighbor_n,
                             verbose = verbose,
                             ...)

  } else if (method == "Seurat"){

    sce <- .integrateSeurat(sce = sce,
                            batch = batch,
                            ndims = ndims,
                            neighbor_n = neighbor_n,
                            verbose = verbose,
                            hvg_ntop = hvg_ntop,
                            ...)

  } else if (method == "LIGER") {

    sce <- .integrateLIGER(sce = sce,
                           batch = batch,
                           ndims = ndims,
                           neighbor_n = neighbor_n,
                           verbose = verbose,
                           ret_nmf = ret_nmf,
                           ...)

  } else if (method == "regression") {

    sce <- .integrateRegression(sce = sce,
                                batch = batch,
                                ndims = ndims,
                                neighbor_n = neighbor_n,
                                verbose = verbose,
                                hvgs = hvgs,
                                ...)

  } else if (method == "STACAS") {
    sce <- .integrateSTACAS(sce = sce,
                            batch = batch,
                            ndims = ndims,
                            neighbor_n = neighbor_n,
                            verbose = verbose,
                            hvg_ntop = hvg_ntop,
                            ...)
  }

  else if (method == "scMerge2") {
    sce <- .integrateScMerge(sce = sce,
                            batch = batch,
                            ndims = ndims,
                            neighbor_n = neighbor_n,
                            verbose = verbose,
                            parallel_param = parallel_param,
                            hvgs = hvgs,
                            ...)
  }

  return(sce)
}

#' @importFrom SummarizedExperiment colData rowData
#' @importFrom uwot umap
#' @importFrom S4Vectors metadata
#' @importFrom batchelor fastMNN

.integrateMNN <- function(sce, batch, hvgs, ndims, neighbor_n, parallel_param, verbose, ...) {

  ## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::integrateSCE() / method = \"MNN\"} - ")

  if (is.null(hvgs)) {
    if (!is.null(metadata(sce)$hvgs)) {
      hvgs <- metadata(sce)$hvgs
    } else stop(ep, "Highly variable genes not supplied and were not previously calculated")
  }

  sce_corr <- fastMNN(sce,
                      batch = colData(sce)[,batch],
                      d = ndims,
                      subset.row = hvgs,
                      BPPARAM = parallel_param,
                      ...)

  reducedDim(sce, "PCA_MNN") <- reducedDim(sce_corr, "corrected")

  # UMAP
  if (verbose) message(.bluem("[INT/fastMNN]"), " Running UMAP on MNN-corrected space")

  reducedDim(sce, "UMAP_MNN") <- umap(reducedDim(sce, "PCA_MNN")[,seq_len(ndims)],
                                      n_neighbors = neighbor_n,
                                      min_dist = 0.7)
  sce
}

#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom uwot umap

.integrateHarmony <- function(sce, batch, dr = "PCA", ndims, neighbor_n, verbose, ...){

  ep <- .redm("{cellula::integrateSCE() / method = \"Harmony\"} - ")

if (!requireNamespace("harmony", quietly = TRUE)) 
   stop(paste0(ep, "the `harmony` package must be installed first.\n
                  Run `BiocManager::install(\"harmony\") to use this function."))
  if(!dr %in% reducedDimNames(sce))
    stop(paste0(ep, "there is no reducedDim slot named ", dr, " in the object."))

  if (verbose) message(.bluem("[INT/Harmony]"), " Correcting batch effect using Harmony")

  harmony_corr <- harmony::RunHarmony(reducedDim(sce, dr)[,seq_len(ndims)],
                                      meta_data = colData(sce)[,batch],
                                      ...)

  reducedDim(sce, "PCA_Harmony") <- harmony_corr

  if (verbose) message(.bluem("[INT/Harmony]"), " Running UMAP on Harmony-corrected space")

  reducedDim(sce, "UMAP_Harmony") <- umap(reducedDim(sce, "PCA_Harmony")[,seq_len(ndims)],
                                          n_neighbors = neighbor_n,
                                          min_dist = 0.7)

  sce
}

#' @importFrom SingleCellExperiment reducedDim counts logcounts
#' @importFrom SummarizedExperiment colData
#' @importFrom uwot umap

.integrateSeurat <- function(sce, batch, ndims, neighbor_n, verbose, hvg_ntop, ...) {

  ep <- .redm("{cellula::integrateSCE() / method = \"Seurat\"} - ")

  if (!requireNamespace("Seurat", quietly = TRUE))
    stop(paste0(ep, "the `Seurat` package must be installed first.\n
                Run `BiocManager::install(\"Seurat\") to use this function."))

  if (verbose) message(.bluem("[INT/Seurat]"), " Converting to Seurat object.")
  old_colnames <- colnames(sce)

  colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))

  nf <- data.frame(old_colnames, row.names = colnames(sce))

  seu <- SeuratObject::CreateSeuratObject(counts = counts(sce),
                                          meta.data = as.data.frame(colData(sce)))
  lc <- logcounts(sce)
  rownames(lc) <- rownames(seu)
  seu <- SeuratObject::SetAssayData(object = seu, layer = "data", new.data = lc)
  seu[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = reducedDim(sce, "PCA"), key = "PC_")

  batches <-  unique(as.character((seu[[batch]][[batch]])))
  seulist <- lapply(batches, function(x) seu[,seu[[batch]] == x])
  names(seulist) <- batches

  if (verbose) message(.bluem("[INT/Seurat]"), " Normalization and HVG selection.")

  seulist <- lapply(seulist, function(x) {
    x <- Seurat::NormalizeData(x, verbose = verbose)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                      nfeatures = hvg_ntop,
                                      verbose = verbose)
  })

  hvgs = Seurat::SelectIntegrationFeatures(seulist, nfeatures = hvg_ntop)

  if (verbose) message(.bluem("[INT/Seurat]"), " Finding anchors.")

  anchors <- Seurat::FindIntegrationAnchors(object.list = seulist,
                                            dims = seq_len(ndims),
                                            verbose = verbose,
                                            anchor.features = hvgs,
                                            ...)

  if (verbose) message(.bluem("[INT/Seurat]"), "Integration.")

  kweight <- min(100, unlist(lapply(seulist, function(x) floor(0.5*(ncol(x))))))

  if (verbose) message(.bluem("[INT/Seurat]"), " Using k.weight = ", kweight, ".")

  seu_int <- Seurat::IntegrateData(anchorset = anchors, dims = seq_len(ndims),
                                   k.weight = kweight)

  if (verbose) message(.bluem("[INT/Seurat]"), " Running PCA on integrated object.")
  seu_int <- Seurat::ScaleData(seu_int,
                               assay = "integrated",
                               verbose = verbose)

  seu_int <- Seurat::RunPCA(seu_int,
                            npcs = ndims,
                            assay = "integrated",
                            verbose = verbose,
                            reduction.name = "spca")

  seu_int <- seu_int[,rownames(nf)]

  if (verbose) message(.bluem("[INT/Seurat]"), " Transferring to SCE object.")

  reducedDim(sce, "PCA_Seurat") <- SeuratObject::Embeddings(seu_int, reduction = "spca")[colnames(sce),]
  rm(seu)
  rm(seu_int)

  if (verbose) message(.bluem("[INT/Seurat]"), " Running UMAP on Seurat-corrected space.")

  reducedDim(sce, "UMAP_Seurat") <- umap(reducedDim(sce, "PCA_Seurat")[,seq_len(ndims)],
                                         n_neighbors = neighbor_n,
                                         min_dist = 0.7)
  colnames(sce) <- nf[colnames(sce), 1]

  sce
}

#' @importFrom SingleCellExperiment reducedDim counts
#' @importFrom SummarizedExperiment colData
#' @importFrom uwot umap

.integrateLIGER <- function(sce, batch, ndims, neighbor_n, ret_nmf, verbose, ...) {

  ep <- .redm("{cellula::integrateSCE() / method = \"LIGER\"} - ")

  if (!requireNamespace("rliger", quietly = TRUE))
    stop(paste0(ep, "the `rliger` and `RcppPlanc` packages must be installed first.\n
                Run `BiocManager::install(c(\"rliger\", \"welch-lab/RcppPlanc\")
                to use this function."))

  old_colnames <- colnames(sce)

  if (any(duplicated(colnames(sce)))) {
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
  }

  nf <- data.frame(old_colnames, row.names = colnames(sce))

  if (verbose) message(.bluem("[INT/LIGER]"), " Creating LIGER object.")

  batches <- unique(colData(sce)[,batch])
  countlist <- lapply(batches, function(x) counts(sce[,colData(sce)[,batch] == x]))
  names(countlist) <- batches
  l <- rliger::createLiger(countlist)
  rm(countlist)

  if (verbose) message(.bluem("[INT/LIGER]"), " Preprocessing.")
  l <- rliger::normalize(l)
  l <- rliger::selectGenes(l)
  l <- rliger::scaleNotCenter(l)

  if (verbose) message(.bluem("[INT/LIGER]"), " Factorization.")
	l <- rliger::runINMF(l, k = neighbor_n,	verbose = verbose, ...)

  if (verbose) message(.bluem("[INT/LIGER]"), " Quantile normalization.")
	l <- rliger::quantileNorm(l, verbose = verbose)		

  if (verbose) message(.bluem("[INT/LIGER]"), " Transferring to SCE object.")
  Hmat = do.call(cbind, lapply(names(l@datasets), function(x) {
  	h = l@datasets[[x]]@H
  	colnames(h) = gsub(paste0("^", x, "_"), "", colnames(h)) 
  	h}))[,colnames(sce)]

  reducedDim(sce, type = "LIGER_iNMF") = t(Hmat)

  for(i in names(l@datasets)) {
	rownames(l@H.norm) = gsub(paste0("^", i, "_"), "",  rownames(l@H.norm)) 
  }
  reducedDim(sce, type = "LIGER_iNMF_NORM") = l@H.norm[colnames(sce),]

  if(ret_nmf) {
	Vmats = lapply(l@datasets, function(x) x@V)
	sce@metadata$LIGER_NMF_matrices = list("H" = Hmat,
                                           "V" = Vmats,
                                            "W" = l@W)
  }

  if (verbose) message(.bluem("[INT/LIGER]"), " Running UMAP on LIGER factorization.")
  reducedDim(sce, "UMAP_LIGER") <- umap(reducedDim(sce, "LIGER_iNMF_NORM")[,seq_len(min(ncol(reducedDim(sce, "LIGER_iNMF_NORM")), ndims))],
                                        n_neighbors = neighbor_n,
                                        min_dist = 0.7)

  colnames(sce) <- nf[colnames(sce), 1]

  sce
}

#' @importFrom SingleCellExperiment reducedDim counts
#' @importFrom SummarizedExperiment colData
#' @importFrom uwot umap
#' @importFrom batchelor regressBatches

.integrateRegression <- function(sce, batch, ndims, hvgs, neighbor_n,
                                 parallel_param, verbose, ...){

  ep <- .redm("{cellula::integrateSCE() / method = \"regression\"} - ")

  if (is.null(hvgs)) {
    if (!is.null(metadata(sce)$hvgs)) {
      hvgs <- metadata(sce)$hvgs
    } else stop(ep, "Highly variable genes not supplied and not previously calculated.")
  }

  if (verbose) message(.bluem("[INT/regression]"), " Correcting batch effect using regression.")

  sce_corr <- regressBatches(sce,
                             batch = colData(sce)[,batch],
                             subset.row = hvgs,
                             d = ndims,
                             BPPARAM = parallel_param,
                             ...)

  reducedDim(sce, "PCA_regression") <- reducedDim(sce_corr, "corrected")

  if (verbose) message(.bluem("[INT/regression]"), " Running UMAP on regression-corrected space.")

  reducedDim(sce, "UMAP_regression") <- umap(reducedDim(sce, "LIGER_NORM")[,seq_len(min(ncol(reducedDim(sce, "PCA_regression")), ndims))],
                                             n_neighbors = neighbor_n,
                                             min_dist = 0.7)
  sce
}

#' @importFrom SummarizedExperiment colData rowData assay assay<-
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom uwot umap
#' @importFrom S4Vectors metadata
#' @importFrom scater runPCA

.integrateScMerge <- function(sce, batch, ndims, hvgs, neighbor_n, parallel_param,
                              verbose, ...) {

  ep <- .redm("{cellula::integrateSCE() / method = \"scMerge2\"} - ")

  if (!requireNamespace("scMerge", quietly = TRUE))
    stop(paste0(ep, "the `scMerge` package must be installed first.\n
                Run `BiocManager::install(\"scMerge\") to use this function."))

  if (verbose) message(.bluem("[INT/scMerge]"), " Correcting batch effect using scMerge2.")

  if (is.null(hvgs)) {
    if (!is.null(metadata(sce)$hvgs)) {
      hvgs <- metadata(sce)$hvgs
    } else stop(ep, "Highly variable genes not supplied and not previously calculated.")
  }

  sce_corr <- scMerge::scMerge2(assay(sce, "logcounts"),
                                batch = colData(sce)[,batch],
                                chosen.hvg = hvgs,
                                verbose = verbose,
                                ...)

  assay(sce, "scMerge2_corrected") <- sce_corr$newY

  if (verbose) message(.bluem("[INT/scMerge2]"), " Running PCA on corrected matrix.")

  sce <- runPCA(sce, exprs_values = "scMerge2_corrected", name = "PCA_scMerge2")

  if (verbose) message(.bluem("[INT/scMerge2]"), " Running UMAP on scMerge2-corrected space.")

  reducedDim(sce, "UMAP_scMerge2") <- umap(reducedDim(sce, "PCA_scMerge2")[,seq_len(min(ncol(reducedDim(sce, "PCA_scMerge2")), ndims))],
                                           n_neighbors = neighbor_n,
                                           min_dist = 0.7)
  sce
}

#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom SingleCellExperiment reducedDim counts logcounts
#' @importFrom uwot umap

.integrateSTACAS <- function(sce, batch, hvg_ntop, ndims, neighbor_n,
                             verbose, ...) {

  ep <- .redm("{cellula::integrateSCE() / method = \"STACAS\"} - ")

  if (!requireNamespace("STACAS", quietly = TRUE))
    stop(paste0(ep, "the `STACAS` package must be installed first.\n
                Run `BiocManager::install(\"carmonalab/STACAS\") to use this function."))

  if (verbose) message(.bluem("[INT/STACAS]"), " Converting to Seurat object.")
  old_colnames <- colnames(sce)

  colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))

  nf <- data.frame(old_colnames, row.names = colnames(sce))

  seu <- SeuratObject::CreateSeuratObject(counts = counts(sce),
                                          meta.data = as.data.frame(colData(sce)))
  batches <-  unique(as.character((seu[[batch]][[batch]])))
  seulist <- lapply(batches, function(x) seu[,seu[[batch]] == x])
  names(seulist) <- batches

  if (verbose) message(.bluem("[INT/STACAS]"), " Normalization.")

  seulist <- lapply(seulist, function(x) {
    x <- Seurat::NormalizeData(x, verbose = verbose)
  })

  hvgs <- Seurat::SelectIntegrationFeatures(seulist, nfeatures = hvg_ntop)

  if(verbose) message(.bluem("[INT/STACAS]"), " Integration.")

  seu_int <- STACAS::Run.STACAS(seulist, anchor.features = hvgs, ...)

  seu_int <- seu_int[,rownames(nf)]

  if (verbose) message(.bluem("[INT/STACAS]"), "Transferring to SCE object.")

  reducedDim(sce, "PCA_STACAS") <- SeuratObject::Embeddings(seu_int, reduction = "pca")[colnames(sce),]
  rm(seu)
  rm(seu_int)

  if (verbose) message(.bluem("[INT/STACAS]"), "Running UMAP on Seurat-corrected space.")

  reducedDim(sce, "UMAP_STACAS") <- umap(reducedDim(sce, "PCA_STACAS")[,seq_len(ndims)],
                                         n_neighbors = neighbor_n,
                                         min_dist = 0.7)
  colnames(sce) <- nf[colnames(sce), 1]

  sce
}


#' Integration diagnostics
#'
#' Diagnose the integration or batch effect correction with different metrics
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name(s) of the column(s) in \code{colData(sce)} with batch labels.
#' 	   It can be more than one column.
#' @param method character, the integration diagnosis method. One of \code{"iLISI"}, \code{"Gini"},
#'     \code{"Jaccard"}, or \code{"HDB"}. Default is \code{"iLISI"}.
#' @param dr character, the name(s) of the dimensional reduction slot(s) to be tested for integration. 
#'     It can be more than one slot. 
#' @param ndims numeric, the number of dimensions to use for integration. Default
#'     is 20. 
#' @param verbose logical, display messages on progress? Default is FALSE.
#' @param ... extra arguments passed to the main diagnosis functions used by each method (see Details).
#'
#' @return depending on the diagnosis methods of choice, a list, a matrix, or a data frame (see Details).
#'
#' @details
#' This function allows to use several dataset integration/batch effect
#' diagnostic methods for single cell datasets, allowing the user to supply
#' any number of parameters using \code{...}.
#' 
#' Depending on the method of choice, the \code{...} placeholder and the output will vary. 
#' 
#' The methods vary as follows:
#' 
#'\itemize{
#'  \item{"iLISI:"}{ integration Local Inverse Simpson's Index from \code{lisi}. 
#'     Provides a per-cell value of "good mixing" within a neighbourhood across one or more
#'     batch variables. It is a numeric index bound between 1 (bad mixing) and the total number
#' 	   of batches for each variable (perfect mixing)}
#' 
#'  \item{"Gini:"}{ Gini index for inequality. Provides a coarser value of "good mixing"
#'     by over-clustering the space using k-means clustering and calculating the Gini index for
#'     the batch variables. Each k-means cluster will be assigned a value between 0 (perfect equality,
#' 	   i.e. perfect mixing) and 1 (perfect inequality, i.e. bad/no mixing).} 
#' 
#'  \item{"Jaccard:"}{ Mean-maximum Jaccard index to compare unbiased clustering and 
#'     ground truth labels. Provides a coarse value of agreement between various clustering resolutions
#' 	   of a space and a previously known cell type annotation, such as one assigned by \code{assignIdentities()}.
#'     
#'     The value is calculated by clustering the data at various sensible resolutions and comparing the
#'     clustering results with ground truth labels, taking the maximum Jaccard index as a measure of 
#' 	   the best possible overlap of a label with any cluster. 
#' 
#'     The mean value of these maximum indices is taken as a measure of the overall 
#'     agreement between the two clusterings. Bound between 0 (no overlap between two 
#'     clusterings) and 1 (perfect overlap between two clusterings). 
#' 	   This method provides some level of information on the preservation of biological variability
#'     after integration, since over-integrated samples should show poor agreement with ground
#'     truth labels by pushing cells with different labels together.}
#' 
#'  \item{"HDB:"}{ Hausdorff Distance of Batches. Provides an asymmetric and statistically
#'     tested pairwise measure of overlap between batches in any space. 
#' 
#'    See the \code{{HDB}} package for more information and arguments to the \code{HDB()} function.}
#'}
#' 
#'  The argument placeholder \code{...} can be used to supply different extra arguments depending on 
#'  the method:
#' \itemize{
#'    \item{"iLISI:"}{ Supply arguments to the
#'     \code{compute_lisi()} function from the \code{{lisi}} package, e.g. \code{perplexity} or
#'     \code{nn_eps}.}
#'    \item{"Gini":}{ Supply the number of k-means clusters using the \code{centers} argument.}
#'    \item{"Jaccard:"}{ Supply the colData column with ground truth labels (\code{labels}), the
#'     numeric vector of clustering resolutions (\code{k}, default = 0.1, 0.3, 0.5, 0.7, 0.9), 
#'     and the number of neighbors for SNN graph construction (\code{neighbors}, default = 30)}
#'    \item{"HDB:"}{ Supply arguments to the \code{HDB} function from the \code{{HDB}} package.}
#' }
#' 
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom SummarizedExperiment colData
#' 
#' @export


diagnoseIntegration <- function(sce,
								batch,
								dr,
								ndims = 20,
								method = "iLISI",
								verbose,
								...){

## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::diagnoseIntegration()} - ")

  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if (!(method %in% c("iLISI", "Gini", "HDB", "Jaccard")))
    stop(ep, "method not recognized - must be one of \"iLISI\", \"Gini\",
         \"HDB\", or \"Jaccard\"")
  if (!(any(batch %in% colnames(colData(sce)))))
    stop(ep,"batch column(s) not found in the colData of the object")
  if (!(any(dr %in% reducedDimNames(sce))))
    stop(ep, "dr slot(s) not found in the object")

  if (method == "iLISI") {
	diagnosis = .diagnoseInt_iLISI(sce = sce, batch = batch, dr = dr, ndims = ndims, verbose = verbose, ...)
  } 
  else if (method == "Gini") {
	diagnosis = .diagnoseInt_Gini(sce = sce, batch = batch, dr = dr, ndims = ndims, verbose = verbose, ...)
  }
  else if (method == "Jaccard") {
	diagnosis = .diagnoseInt_Jaccard(sce = sce, batch = batch, dr = dr, ndims = ndims, verbose = verbose, ...)
  }
  else if (method == "HDB") {
	diagnosis = .diagnoseInt_HDB(sce = sce, batch = batch, dr = dr, ndims = ndims, verbose = verbose, ...)
  }
  return(diagnosis)
}


#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames

.diagnoseInt_iLISI <- function(sce, batch, dr, ndims, verbose, ...) {
	
## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::diagnoseIntegration() / method = \"iLISI\"} - ")

 if (!requireNamespace("lisi", quietly = TRUE)) 
   		stop(paste0(ep, "the `lisi` package must be installed first.\n
                  Run `BiocManager::install(\"immunogenomics/lisi\") to use this function."))

	lisil = lapply(dr, 
                function(r) {
                  if(verbose) message(paste0(.bluem("[INT/Diagnosis - iLISI]"), " Calculating iLISI for reduction ", r))
                  lisi::compute_lisi(reducedDim(sce, r)[,seq_len(ndims)], 
				  meta_data = colData(sce), 
				  label_colnames = batch,
				  ...)
                })
	if(length(dr) > 1) {
		names(lisil) = dr
		out = do.call(cbind, lisil)
		rownames(out) = colnames(sce)
	} else {
		out = unlist(lisil)
		names(out) = colnames(sce)
	}	
	return(out)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom stats kmeans
#' @importFrom utils combn

.diagnoseInt_Gini <- function(sce, batch, dr, ndims, verbose, centers = 200) {

## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::diagnoseIntegration() / method = \"Gini\"} - ")

  ginil = lapply(dr,
    function(r) {
    	if(verbose) message(paste0(.bluem("[INT/Diagnosis - Gini]"), " Calculating Gini index for reduction ", r))
         sce_overcl = kmeans(reducedDim(sce, r)[,seq_len(ndims)], centers = centers)$cluster
         index_tab = table(sce_overcl, colData(sce)[,batch])
         unlist(lapply(seq_len(nrow(index_tab)), function(x) {
           diffs = abs(apply(t(combn(index_tab[x,], 2)), 1, diff))
           gini_int = ((sum(diffs)/length(diffs))/2)/mean(index_tab[x,])
		   gini_int
         }))
       })
if(length(dr) > 1) {
		names(ginil) = dr
		out = do.call(cbind, ginil)
	} else {
		out = unlist(ginil)
	}	
	return(out)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDimNames

.diagnoseInt_Jaccard <- function(sce, batch, dr, ndims, verbose, labels, neighbors = 30, k = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
	
	## Sanity checks
	# Error prefix
	ep <- .redm("{cellula::diagnoseIntegration() / method = \"Jaccard\"} - ")

	if (!(labels %in% colnames(colData(sce))))
		stop(ep, "labels must be a column of the colData")
		if(length(k) < 1)
		stop(ep, "specify at least 1 value of k")

	jacl = lapply(dr,
		function(r) {
    	if(verbose) message(paste0(.bluem("[INT/Diagnosis - Jaccard]"), " Calculating mean-max Jaccard index for reduction ", r))
			
			sce_overcl = colData(makeGraphsAndClusters(sce, neighbors = neighbors, k = k, dr = r))[,paste0("SNN_", k)]
			
			jac = lapply(sce_overcl, function(x) {
				pairs = as.data.frame(expand.grid(unique(x), unique(colData(sce)[,labels])))
				pairs$intersect = apply(pairs, 1, function(y)
				length(intersect(which(x == y[1]), which(colData(sce)[,labels] == y[2]))))
				pairs$union = apply(pairs, 1, function(y) {
					length(union(which(x == y[1]), which(colData(sce)[,labels] == y[2])))
				})
				pairs$jaccard = pairs$intersect/pairs$union
				jacmat = matrix(pairs$jaccard, ncol = length(unique(colData(sce)[,labels])))
				mmj = mean(colMaxs(jacmat))
				return(list("jacmat" = jacmat, "meanmax_jaccard" = mmj))
			}
		)
		list(clusters = sce_overcl,jaccard = jac)
		}
	)
	if(length(dr) > 1) {
		names(jacl) = dr
		out = do.call(rbind, lapply(jacl, function(x) 
   			unlist(lapply(x$jaccard, function(y) y$meanmax_jaccard))))
	} else {
		out = unlist(jacl)
	}
	return(out)
}

#' @importFrom BiocParallel SerialParam

.diagnoseInt_HDB <- function(sce, batch, dr, ndims, verbose, ...) {

	## Sanity checks
	# Error prefix
	ep <- .redm("{cellula::diagnoseIntegration() / method = \"HDB\"} - ")

	if (!requireNamespace("HDB", quietly = TRUE)) 
			stop(paste0(ep, "the `HDB` package must be installed first.\n
					Run `BiocManager::install(\"gdagstn/HDB\") to use this function."))
	
	if (!requireNamespace("BiocNeighbors", quietly = TRUE)) 
			stop(paste0(ep, "the `HDB` package must be installed first.\n
					Run `BiocManager::install(\"BiocNeighbors\") to use this function."))
	hdl = lapply(dr, function(r) {
	    	if(verbose) message(paste0(.bluem("[INT/Diagnosis - HDB]"), " Calculating HDB for reduction ", r))

		HDB::HDB(sce = sce, dimred = r, group = batch, distance = "Euclidean", q = 3, dims = ndims, 
			 	props = c(0.2, 0.5, 0.8), samples = 50, BPPARAM = SerialParam(), doplot = FALSE, verbose = verbose, ...)
	})	
	names(hdl) = dr
	return(hdl)
}
