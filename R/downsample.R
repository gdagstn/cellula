#' Downsample counts
#'
#' Downsamples a count matrix to a predetermined total count
#'
#' @param mat a matrix where rows are genes/features and columns are samples
#' @param target numeric, the total number of counts per sample to be reached by
#'     downsampling
#' @param chunksize numeric, the size (in number of samples) to be processed at
#'     a time
#' @param BPPARAM a BiocParallel BPPARAM object. Default is SerialParam(), meaning
#'     no parallelization will be used.
#' @param verbose logical, should progress be displayed? Default is TRUE.
#'
#' @return a sparse matrix of class "dgCMatrix"
#'
#' @author Giuseppe D'Agostino, inspired by Scott R. Tyler's \code{downsample} package
#'
#' @importFrom methods as
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom Matrix sparseVector

downsampleCounts <- function(mat,
                             target = NULL,
                             chunksize = 1000,
                             BPPARAM = SerialParam(progressbar = verbose),
                             verbose = TRUE) {

  all_genes = rownames(mat)

  if(is.null(target)) target = min(colSums(mat))

  if(!is.null(target)) {
    if(all(colSums(mat) < target)) {
      stop("No cells left above target!")
    }
    if(any(colSums(mat) < target)) {
      if(verbose) message("Removing cells with total count below target")
      mat = mat[,colSums(mat) >= target]
    }
  }

  chunks = split(seq_len(ncol(mat)),
                 ceiling(seq_along(seq_len(ncol(mat)))/chunksize))

  ds = bplapply(chunks, function(n) {
    ch = mat[,n]
    dsch = apply(ch, 2, function(x) {
      nz_genes = rownames(mat)[which(x > 0)]
      nz_genes = rep(nz_genes, x[nz_genes])

      sim = sample(nz_genes,
                   size = target,
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

  ds_final = do.call(cbind, ds)
  ds_final
}

