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
#' @importFrom crayon blue
#' @importFrom uwot umap
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom BiocParallel SerialParam
#' @importFrom scater runPCA
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export
doNormAndReduce <- function(sce, batch = NULL, name = NULL,
                            ndims = 20,
                            hvg_ntop = 2000,
                            verbose = TRUE,
                            parallel_param = SerialParam()) {
  
  #Sanity checks
  # Error prefix
  ep = "{cellula::doNormAndReduce()} - "
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "must provide a SingleCellExperiment object"))
  if(!is.null(batch)) {
    if(!batch %in% colnames(colData(sce)))
      stop(paste0(ep, "batch column not found in the colData of the object"))
  }
  if(hvg_ntop > nrow(sce))
    stop(paste0(ep, "hvg_ntop cannot be higher than the number of features (nrow) in the object"))
  
  if(is.null(metadata(sce)$cellula_log)) {
    clog = .initParams()
  } else {
    clog = metadata(sce)$cellula_log
  }
  
  if(verbose) cat(blue("[NORM]"), "Calculating size factors and normalizing. \n")

  # Size factors
  if(verbose) cat(blue("[NORM]"),"   Preclustering. \n")

  if(!is.null(batch)) {
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

  if(verbose) cat(blue("[NORM]"),"   Calculating pooled factors. \n")
  sce <- computeSumFactors(sce,
                           clusters = sce_cl,
                           BPPARAM = parallel_param)

  # Normalization
  if(verbose) cat(blue("[NORM]"),"   Log-normalization. \n")
  if(!is.null(batch)) {
    sce <- multiBatchNorm(sce,
                          batch = colData(sce)[,batch])
    clog$norm_reduce$norm_strategy = "multiBatchNorm"
  } else {
    sce <- logNormCounts(sce)
    clog$norm_reduce$norm_strategy = "logNormCounts"
  }
  
  if(!is.null(name)) {
    if(verbose) cat("Saving temporary file. \n")
    saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))
  }
  # HVGs

  if(verbose) cat(blue("[DR]"), "Selecting HVGs. \n")

  if(!is.null(batch)) {
    vargenes <- modelGeneVar(sce,
                        block = colData(sce)[,batch])
  } else {
    vargenes <- modelGeneVar(sce)
  }

  hvgs <- getTopHVGs(vargenes, n = hvg_ntop)
  metadata(sce)$hvgs <- hvgs
  
  clog$norm_reduce$hvg_ntop = hvg_ntop
  
  # PCA

  if(verbose) cat(blue("[DR]"), "Running PCA. \n")

  sce <- runPCA(sce,
                subset_row = hvgs,
                exprs_values = "logcounts",
                ncomponents	= ndims)#,
  #BPPARAM = parallel_param)

  if(!is.null(name)) {
    if(verbose) cat("Saving temporary file. \n")
    saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))
  }

  # UMAP
  if(verbose) cat(blue("[DR]"), "Running UMAP on uncorrected PCA \n")

  neighbor_n <- floor(sqrt(ncol(sce)))

  reducedDim(sce, "UMAP") <- umap(reducedDim(sce, "PCA")[,seq_len(ndims)],
                                  n_neighbors = neighbor_n,
                                  min_dist = 0.7)
  
  clog$norm_reduce$umap_min_dist = 0.7
  clog$norm_reduce$umap_n_neighbors = neighbor_n
  clog$norm_reduce$umap_other = formals(umap)[!formalArgs(umap) %in% c("n_neighbors", "min_dist", "X")]
  clog$norm_reduce$parallel_param = parallel_param
  
  metadata(sce)$cellula_log = clog
    
  return(sce)
}
