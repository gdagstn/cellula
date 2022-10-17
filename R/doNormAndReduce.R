
#' SCE Normalization and dimensionality reduction sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels
#' @param name character, name of the folder and file where results are stored
#' @param ndims numeric, number of dimensions (PCs) to be used for UMAP construction.
#'     Default is 20.
#' @param hvg_ntop numeric, number of top HVGs to retain
#' @param verbose logical, display messages on progress? Default is FALSE.
#' @param parallel_param a BiocParallel object specifying the parallelization backend
#'     to be used in some steps of the pipeline.
#'
#' @return a SingleCellExperiment object containing normalized values and a non-integrated
#'    PCA and UMAP dimensionality reductions.
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom scran quickCluster computeSumFactors modelGeneVar getTopHVGs
#' @importFrom batchelor multiBatchNorm
#' @importFrom crayon blue
#' @importFrom uwot umap
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors metadata metadata<-

doNormAndReduce <- function(sce, batch, name,
                            ndims = 20,
                            hvg_ntop = 2000,
                            verbose = TRUE,
                            parallel_param = SerialParam()) {

  if(verbose) cat(blue("[NORM]"), "Calculating size factors and normalizing. \n")

  # Size factors
  if(verbose) cat(blue("[NORM]"),"   Preclustering. \n")

  sce_cl <- quickCluster(sce,
                         block = colData(sce)[,batch],
                         min.size = min(100, min(table(colData(sce)[,batch]))),
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
  metadata(sce)$hvgs <- hvgs

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

  return(sce)
}
