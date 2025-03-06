#' SCE empty droplet removal sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels.
#'     Default is NULL meaning no batches will be considered, and data will be
#'     processed as a single batch.
#' @param emptydrops_cutoff either \code{"auto"} (default, barcode rank inflection point)
#'     or a numeric. Cells with total reads below this cutoff are used to calculate
#'     ambient RNA profiles and are removed.
#' @param emptydrops_alpha numeric, the FDR threshold to call an empty barcode.
#' @param verbose logical, display messages on progress? Default is \code{FALSE}.
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Default is \code{SerialParam()} 
#'     meaning no parallelization will be used.
#'
#' @return a \code{SingleCellExperiment} object without empty droplets as assigned by \code{emptyDrops}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export

doEmptyDrops <- function(sce,
                          batch = NULL,
                          emptydrops_cutoff = "auto",
                          emptydrops_alpha = 0.01,
                          verbose = TRUE,
                          parallel_param = SerialParam()) {

  ep <- .redm("{cellula::doEmptyDrops()} - ")
  
	dependencies = data.frame("package" = c("DropletUtils"),
                              "repo" = c("BioC"))
    if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))

  
  if (verbose) message(.bluem("[QC/EMPTY]"),"Running emptyDrops.")

  if (!is.null(batch)) {
    batches <- unique(as.character(colData(sce)[,batch]))
    scelist <- lapply(batches, function(x) sce[,colData(sce)[,batch] == x])
    names(scelist) <- batches
    if (emptydrops_cutoff == "auto") {
      emptydrops_cutoff <- lapply(scelist, function(x) 
        metadata(DropletUtils::barcodeRanks(x))$inflection)
    } else {
      emptydrops_cutoff <- as.list(rep(emptydrops_cutoff, length(scelist)))
    }

    names(emptydrops_cutoff) <- batches

    empty_list <- bplapply(names(scelist), function(x) {
                            empty_droplets <- DropletUtils::emptyDrops(scelist[[x]],
                                                         lower = emptydrops_cutoff[[x]])
                            keep_droplets <- empty_droplets$FDR <= emptydrops_alpha
                            scelist[[x]]$empty <- factor(ifelse(empty_droplets$FDR <= emptydrops_alpha, "ok", "empty"))
                            scelist[[x]]$empty[which(is.na(scelist[[x]]$empty))] <- "empty"
                            if (verbose) message(table(Sig = keep_droplets, Limited = empty_droplets$Limited))
                            if (verbose) message("Empty cells (", x, "): ", sum(scelist[[x]]$empty == "empty"))
                            return(scelist[[x]]$empty)
                          }, BPPARAM = parallel_param)

    names(empty_list) <- names(scelist)
    sce$empty <- unlist(empty_list)
   } else {
    if (emptydrops_cutoff == "auto") {
      barcode_ranks <- DropletUtils::barcodeRanks(sce)
      emptydrops_cutoff <- metadata(barcode_ranks)$inflection
    }

    empty_droplets <- DropletUtils::emptyDrops(sce, lower = emptydrops_cutoff)
    keep_droplets <- empty_droplets$FDR <= emptydrops_alpha
    sce$empty <- factor(ifelse(empty_droplets$FDR <= emptydrops_alpha, "ok", "empty"))
    sce$empty[which(is.na(sce$empty))] <- "empty"

    if (verbose) print(table(Sig = keep_droplets, Limited = empty_droplets$Limited))
    if (verbose) message("Empty cells: ", sum(sce$empty == "empty"))
  }
  sce <- sce[, which(sce$empty == "ok"), drop = FALSE]
  sce
}
