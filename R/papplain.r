#'#' SCE Integration pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels
#' @param name character, the name of the project/folder where files will be saved.
#'     If NULL, a random name will be generated
#' @param do_qc logical, should QC steps be performed? Default is TRUE
#' @param do_norm logical, should normalization and dimensionality reduction be
#'     performed? Default is TRUE
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
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom BiocParallel SerialParam
#'
#' @export

papplain <- function(sce,
                     batch = NULL,
                     name = NULL,
                     do_qc = TRUE,
                     do_norm = TRUE,
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

 if(!is.null(batch)) {
   if(!batch %in% colnames(colData(sce)))
    stop(paste0("Batch label \"", batch, "\" not found."))
 }

  # Make folder
  if(is.null(name)) {
    name = adjective_animal()
    cat(paste0("No name selected so the randomly assigned name is: ", name, "\n"))
  }

  dir.create(name)
 if(!is.null(batch)) {
   if(!is.factor(colData(sce)[,batch]))
     colData(sce)[,batch] = as.factor(colData(sce)[,batch])
 }

  if(verbose) {
    cat("Working on object", name, "\n")
    ncells = ncol(sce)
    cat("Input cells: ", ncells, "\n")
    if(!is.null(batch)) {
      cat("By batch:\n")
      print(table(colData(sce)[,batch]))
    }
  }

  if(do_qc) {
    sce <- doQC(sce, name = name,
                batch = batch, discard = discard,
                subset_mito = subset_mito,
                subset_ribo = subset_ribo,
                subset_malat1 = subset_malat1,
                detect_doublets = detect_doublets,
                run_emptydrops = run_emptydrops,
                emptydrops_cutoff = emptydrops_cutoff,
                emptydrops_alpha = emptydrops_alpha,
                verbose = verbose,
                save_plots = save_plots,
                parallel_param = parallel_param)

    #if(!is.null(stopat) & stopat == "DBL") return(sce)
    if(verbose) cat("Saving temporary file. \n")

    saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))
  }

  if(do_norm) {

    sce <- doNormAndReduce(sce, batch = batch, name = name,
                           ndims = ndims,
                           hvg_ntop = hvg_ntop,
                           verbose = verbose,
                           parallel_param = parallel_param)

    if(verbose) cat("Saving temporary file. \n")
    saveRDS(sce, file = paste0("./", name, "/", name, "_tempSCE.RDS"))
  }

  if(verbose) cat(blue("[INT]"), "Integration. \n")
  hvgs = metadata(sce)$hvgs
  if(!is.null(batch)) {
    sce = integrateSCE(sce,
                       batch = batch,
                       hvgs = hvgs,
                       hvg_ntop = hvg_ntop,
                       method = integration_method,
                       ndims = ndims,
                       parallel_param = parallel_param,
                       verbose = verbose)
  }

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



