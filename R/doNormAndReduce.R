#' SCE Normalization and dimensionality reduction sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels.
#'     Default is \code{NULL} meaning no batches will be considered, and data will be
#'     processed as a single batch.
#' @param name character, name of the folder and file where results are stored
#' @param ndims numeric, number of dimensions (PCs) to be used for UMAP construction.
#'     Default is 20.
#' @param hvg_ntop numeric, number of top HVGs to retain
#' @param verbose logical, display messages on progress? Default is \code{FALSE}.
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Default is \code{SerialParam()},
#'     meaning no parallelization will be used.
#'
#' @return a \code{SingleCellExperiment} object containing normalized values and a non-integrated
#'    \code{PCA} and \code{UMAP} dimensionality reductions.
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom scran quickCluster computeSumFactors modelGeneVar getTopHVGs
#' @importFrom scuttle logNormCounts
#' @importFrom batchelor multiBatchNorm
#' @importFrom uwot umap
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom BiocParallel SerialParam
#' @importFrom scater runPCA
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom methods formalArgs
#'
#' @export

doNormAndReduce <- function(sce, batch = NULL, name = NULL,
                            ndims = 20,
                            hvg_ntop = 2000,
                            verbose = TRUE,
                            parallel_param = SerialParam()) {
  
  # Checks
  ep <- .redm("{cellula::doNormAndReduce()} - ")
  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "must provide a SingleCellExperiment object")
  if (!is.null(batch)) {
    if (!batch %in% colnames(colData(sce)))
      stop(ep, "batch column not found in the colData of the object")
  }
  if (hvg_ntop > nrow(sce))
    stop(ep, "hvg_ntop cannot be higher than the number of features (nrow) in the object")
  # Start parameter logging - not fully implemented
  if (is.null(metadata(sce)$cellula_log)) {
    clog <- .initLog()
  } else {
    clog <- metadata(sce)$cellula_log
  }
    if (verbose) message(.bluem("[NORM] "), "Calculating size factors and normalizing.")
  # Size factors
  if (verbose) message(.bluem("[NORM] "),"   Preclustering.")
  if (!is.null(batch)) {
    sce_cl <- quickCluster(sce,
                         block = colData(sce)[,batch],
                         min.size = floor(sqrt(min(table(colData(sce)[,batch])))),
                         BPPARAM = parallel_param)
    clog$norm_reduce$precluster_min_size = floor(sqrt(min(table(colData(sce)[,batch]))))
  } else {
    sce_cl <- quickCluster(sce,
                           min.size = floor(sqrt(ncol(sce))),
                           BPPARAM = parallel_param)
    clog$norm_reduce$precluster_min_size = floor(sqrt(ncol(sce)))
  }
  if (verbose) message(.bluem("[NORM] "),"   Calculating pooled factors.")
  sce <- computeSumFactors(sce,
                           clusters = sce_cl,
                           BPPARAM = parallel_param)
  # Normalization
  if (verbose) message(.bluem("[NORM] "),"   Log-normalization.")
  if (!is.null(batch)) {
    sce <- multiBatchNorm(sce,
                          batch = colData(sce)[,batch])
    clog$norm_reduce$norm_strategy = "multiBatchNorm"
  } else {
    sce <- logNormCounts(sce)
    clog$norm_reduce$norm_strategy = "logNormCounts"
  }
  
  # HVGs
  if (verbose) message(.bluem("[DR] "), "Selecting HVGs.")

  if (!is.null(batch)) {
    vargenes <- modelGeneVar(sce,
                        block = colData(sce)[,batch])
  } else {
    vargenes <- modelGeneVar(sce)
  }
  hvgs <- getTopHVGs(vargenes, n = hvg_ntop)
  metadata(sce)$hvgs <- hvgs
  clog$norm_reduce$hvg_ntop = hvg_ntop
  # PCA
  if(verbose) message(.bluem("[DR] "), "Running PCA.")
  sce <- runPCA(sce,
                subset_row = hvgs,
                exprs_values = "logcounts",
                ncomponents	= ndims)#,
  #BPPARAM = parallel_param) # the overhead for parallel PCA seems to be big.
  # UMAP
  if (verbose) message(.bluem("[DR] "), "Running UMAP on uncorrected PCA.")
  neighbor_n <- floor(sqrt(ncol(sce)))
  reducedDim(sce, "UMAP") <- umap(reducedDim(sce, "PCA")[,seq_len(ndims)],
                                  n_neighbors = neighbor_n,
                                  min_dist = 0.7)
  # Parameter logs
  clog$norm_reduce$umap_min_dist <- 0.7
  clog$norm_reduce$umap_n_neighbors <- neighbor_n
  clog$norm_reduce$umap_other <- formals(umap)[!formalArgs(umap) %in% c("n_neighbors", "min_dist", "X")]
  clog$norm_reduce$parallel_param <- parallel_param
  
  metadata(sce)$cellula_log <- clog
    
  sce
}


#' UMAP Parameter Sweep
#'
#' Runs 
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param dr character, name of the \code{reducedDim} slot in \code{sce}
#' @param ndims numeric, number of dimensions (PCs) to be used for UMAP construction.
#'     Default is 20.
#' @param md numeric vector, values for the minimum distance parameter grid
#' @param nn numeric vector, values for the number of neighbors parameter grid
#' @param seed numeric, the random number generator seed
#' @param consistent_seed logical, should the same seed be applied to all runs? 
#'     Default is FALSE
#' @param verbose logical, display messages on progress? Default is \code{FALSE}.
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Default is \code{SerialParam()},
#'     meaning no parallelization will be used.
#'
#' @return a \code{list} containing two elements: a \code{list} of UMAP embedding coordinates and a \code{data.frame} with the parameter grid and its associated seed
#'
#' @importFrom uwot umap
#' @importFrom SingleCellExperiment reducedDim 
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @export

paramSweepUMAP <- function(sce, dr, ndims = 20, md = seq(0.1, 1, length.out = 10), 
                           nn = seq(5, 100, length.out = 20), seed = 11, 
                           consistent_seed = FALSE, verbose = FALSE,
                           parallel_param = SerialParam())
  {

    set.seed(seed)

    if(consistent_seed) 
    {
      rng_seeds = rep(seed, length(md) * length(nn)) 
    } else {
      rng_seeds = sample(size = length(md) * length(nn), x = seq_len(1e10))
    }
    
    combs = as.data.frame(expand.grid(md, nn))
    comps$seed = rng_seeds

    space = reducedDim(sce, dr)[,seq_len(ndims)]

    umaps = bplapply(seq_len(nrow(combs)), function(x) {
      umap(space, min_dist = combs[x,1], n_neighbors = combs[x,2], seed = combs[x,3])
    }, BPPARAM = parallel_param)

    return(list("umaps" = umaps, "parameters" = combs))
}