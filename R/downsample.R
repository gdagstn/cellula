#' Downsample counts
#'
#' Downsamples a count matrix to a predetermined total count
#'
#' @param sce a \code{SingleCellExperiment} object from which the \code{counts}
#'     assay will be taken.
#' @param target numeric, the total number of counts per sample to be reached by
#'     downsampling
#' @param chunksize numeric, the size (in number of samples) to be processed at
#'     a time
#' @param BPPARAM a \code{BiocParallel} \code{BPPARAM} object. 
#'     Default is \code{SerialParam()}, meaning no parallelization will be used.
#' @param verbose logical, should progress be displayed? Default is \code{TRUE}.
#'
#' @return a sparse matrix of class \code{"dgCMatrix"} with sampled counts
#'
#' @author Giuseppe D'Agostino, inspired by Scott R. Tyler's \code{downsample} package
#'
#' @importFrom methods as
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom Matrix sparseVector colSums
#' @importFrom SummarizedExperiment assay
#' 
#' @export

downsampleCounts <- function(sce,
                             target = NULL,
                             chunksize = 1000,
                             BPPARAM = SerialParam(progressbar = verbose),
                             verbose = TRUE) {
  
  # Checks
  ep = .redm("{cellula::downsampleCounts()} - ")
  if (!is(sce, "SingleCellExperiment")) 
    stop(ep, "must provide a SingleCellExperiment object")
  if (!is.null(target) & !is(target, numeric)) 
    stop(ep, "target must be a numeric or NULL")
  
  # Start parameter logging - not fully implemented
  # TO DO
  # --------------------------------------------- #
  
  mat = assay(sce, "counts")
  all_genes = rownames(mat)

  if (is.null(target)) target = min(colSums(mat))
  if (!is.null(target)) {
    if (all(colSums(mat) < target)) {
      stop("No cells left above target!")
    }
    if(any(colSums(mat) < target)) {
      if(verbose) message("Removing cells with total count below target")
      mat = mat[,colSums(mat) > target]
    }
  }

  chunks = split(seq_len(ncol(mat)),
                 ceiling(seq_along(seq_len(ncol(mat)))/chunksize))

  ds = bplapply(chunks, function(n) {
    ch = mat[,n,drop=FALSE]
    dsch = apply(ch, 2, function(x) {
      nz_genes = rownames(mat)[which(x > 0)]
      nz_genes = rep(nz_genes, x[nz_genes])

      sim = sample(nz_genes,
                   size = min(target, length(nz_genes)),
                   replace = FALSE)

      genes_sampled = table(sim)[nz_genes]
      genes_sampled = genes_sampled[!is.na(genes_sampled)]
      x[all_genes] = 0
      x[names(genes_sampled)] = genes_sampled
      names(x) = NULL
      xs = sparseVector(x[x>0],
                        length = length(x),
                        i = which(x > 0))
      xs
    })
    dsc = as(do.call(cbind, lapply(dsch, as.matrix, sparse = TRUE)), "dgCMatrix")
    dimnames(dsc) = dimnames(mat[,n])
    dsc
  },
  BPPARAM = BPPARAM)
  do.call(cbind, ds)
}

#' Downsample cells
#'
#' Downsamples a SingleCellExperiment to a proportion of cells 
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param sample_by character, the column name of \code{colData(sce)} for sampling. 
#'     Cells will be downsampled within each level of \code{sample_by}.
#' @param min numeric, the minimum number of cells to retain within each level of
#'     \code{sample_by}. Strictly positive.
#' @param proportion numeric, the target proportion of cells to downsample within
#'     each level of \code{sample_by}. 
#' @param BPPARAM a \code{BiocParallel} \code{BPPARAM} object. 
#'     Default is \code{SerialParam()}, meaning no parallelization will be used.
#' @param verbose logical, should progress be displayed? Default is \code{TRUE}.
#'
#' @return a \code{SingleCellExperiment} object with fewer cells than the original
#'
#' @author Giuseppe D'Agostino
#'
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom SummarizedExperiment colData
#' 
#' @export

downsampleCells <- function(sce, 
                            sample_by, 
                            proportion = 0.1, 
                            min = 10, 
                            verbose = TRUE, 
                            BPPARAM = SerialParam()) {
  
  ## Sanity checks
  # Error prefix
  ep = .redm("{cellula::downsampleCells()} - ")
  
  if (!is(sce, "SingleCellExperiment")) 
    stop(ep, "must provide a SingleCellExperiment object")
  if(min <= 0) 
    stop(ep, "min must be strictly positive")
  if(!is(min, "numeric")) 
    stop(ep, "min must be a numeric")
  if(!is(proportion, "numeric")) 
    stop(ep, "proportion must be a numeric")
  if(proportion > 1 | proportion < 0) 
    stop(ep, "proportion must be between 0 and 1 excluded")
  if(!sample_by %in% colnames(colData(sce))) 
    stop(ep, sample_by, " is not a column in colData(sce)")

  sampled = bplapply(unique(colData(sce)[,sample_by]), function(x) {
    if(verbose) message(.bluem("[DOWNSAMPLE]"), "Downsampling cells in ", x)
    curr = sce[, colData(sce)[,sample_by] == x]
    if(ncol(curr) < min) {
      warning("Total cell number for ", x, " is less than min. \nWill retain all cells in ", x, ".")
      min = ncol(curr)
    }
    subp = max(c(floor(ncol(curr)*proportion), min))
    keep = colnames(curr)[sample(x = seq_along(colnames(curr)), size = subp)]
  },
  BPPARAM = BPPARAM)
  keep = Reduce(c, sampled)
  sce[,keep]
}

