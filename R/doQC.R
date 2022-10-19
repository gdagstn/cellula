#' SCE QC sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels.
#'     Default is NULL meaning no batches will be considered, and data will be
#'     processed as a single batch.
#' @param name character, the name of the file/folder.
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
#' @param verbose logical, display messages on progress? Default is FALSE.
#' @param save_plots logical, should plots be drawn and saved? Default is TRUE
#' @param parallel_param a BiocParallel object specifying the parallelization backend
#'     to be used in some steps of the pipeline.
#'
#' @return  a `SingleCellExperiment` object with QC including doublet assignment
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom DropletUtils emptyDrops barcodeRanks
#' @importFrom crayon blue
#' @importFrom scuttle isOutlier perCellQCMetrics
#' @importFrom gridExtra grid.arrange
#' @importFrom scater plotColData
#' @importFrom ggplot2 scale_y_log10 ggtitle ggsave
#' @importFrom scDblFinder scDblFinder
#' @importFrom BiocParallel SerialParam
#'
#' @export

doQC <- function(sce,
                 batch = NULL,
                 name,
                 discard = TRUE,
                 subset_mito = TRUE,
                 subset_ribo = TRUE,
                 subset_malat1 = TRUE,
                 detect_doublets = TRUE,
                 run_emptydrops = FALSE,
                 emptydrops_cutoff = "auto",
                 emptydrops_alpha = 0.01,
                 verbose = TRUE,
                 save_plots = TRUE,
                 parallel_param = SerialParam()){

  if(run_emptydrops){

    sce <- doEmptyDrops(sce, batch = batch,
                        emptydrops_cutoff = emptydrops_cutoff,
                        emptydrops_alpha = emptydrops_alpha,
                        verbose = verbose,
                        parallel_param = parallel_param)
  }

  if(verbose) cat(blue("[QC]"),"Calculating QC metrics. \n")

  # Generate QC
  if(subset_mito){
    if(sum(grepl("^MT-", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      cat("   No MT genes found.\n")
      mito = FALSE
    } else mito = rownames(sce)[grepl("^MT-", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else mito = FALSE

  if(subset_malat1){
    if(sum(grepl("^MALAT1", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      cat("   No MALAT1 gene found.\n")
      Malat1 = FALSE
    } else Malat1 = rownames(sce)[grepl("^MALAT1", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else Malat1 = FALSE

  if(subset_ribo){
    if(sum(grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      cat("   No ribo genes found.\n")
      Ribo = FALSE
    } else Ribo = rownames(sce)[grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else Ribo = FALSE

  subset_list = list(mito = mito, Malat1 = Malat1, Ribo = Ribo)
  subset_list = subset_list[!is.na(subset_list)]

  if(length(subset_list) == 0) subset_list = NULL

  if(!is.null(batch)) qcbatch = colData(sce)[,batch] else qcbatch = NULL

    sce_fqc <- perCellQCMetrics(sce,
                                subsets = subset_list,
                                BPPARAM = parallel_param)

    low.lib <- isOutlier(log10(sce_fqc$sum),
                         batch = qcbatch,
                         type = "lower",
                         nmads=3)

    low.genes <- isOutlier(log10(sce_fqc$detected),
                           batch = qcbatch,
                           type = "lower",
                           nmads=3)

    if(!all(is.na(mito))) {
      high.mt <- isOutlier(sce_fqc$subsets_mito_percent,
                           batch = qcbatch,
                           type = "higher",
                           nmads = 3)
    } else high.mt <- FALSE

    if(!all(is.na(Malat1))) {
      high.malat1 <- isOutlier(sce_fqc$subsets_Malat1_percent,
                               batch = qcbatch,
                               type = "higher",
                               nmads = 3)
    } else high.malat1 <- FALSE

    if(!all(is.na(Ribo))) {
      high.ribo <- isOutlier(sce_fqc$subsets_Ribo_percent,
                             batch = qcbatch,
                             type = "higher",
                             nmads = 3)
    } else high.ribo <- FALSE

    data.frame(LowLib=sum(low.lib),
               LowNgenes=sum(low.genes),
               HighMT = ifelse(all(mito == FALSE), FALSE, sum(high.mt)),
               HighMalat1 = ifelse(all(Malat1 == FALSE), FALSE, sum(high.malat1)),
               HighRibo = ifelse(all(Ribo == FALSE), FALSE, sum(high.ribo)))

    sce_fqc$discard <- low.lib | low.genes | high.mt | high.malat1 | high.ribo

    colData(sce) <- cbind(colData(sce), sce_fqc)

    if(save_plots){

      if(verbose) cat(blue("[QC]"),"   Saving QC plots. \n")

      p1 <- grid.arrange(
        plotColData(sce, y="sum", x = batch, colour_by="discard") +
          scale_y_log10() + ggtitle("Total count"),
        plotColData(sce, y="detected", x = batch, colour_by="discard") +
          scale_y_log10() + ggtitle("Detected features"),
        plotColData(sce, y="subsets_Malat1_percent", x = batch,
                    colour_by="discard") + ggtitle("Malat1 percent"),
        plotColData(sce, y="subsets_Ribo_percent", x = batch,
                    colour_by="discard") + ggtitle("Ribo percent"),
        plotColData(sce, y="subsets_mito_percent", x = batch,
                    colour_by="discard") + ggtitle("Mito percent"),
        ncol=2
      )
      ggsave(p1, filename = paste0("./", name, "/QC_plot.png"), width =  6, height = 6, device = "png")
    }

    sce$discard[is.na(sce$discard)] = TRUE

    if(discard) sce <- sce[, !sce$discard]

  # Doublet finding
  if(detect_doublets){
    if(verbose) cat(blue("[QC/DBL]"), "Finding doublets. \n")
    if(!is.null(batch)) samples = colData(sce)[,batch] else samples = NULL
    sce <- scDblFinder(sce, verbose = verbose,
                       samples = colData(sce)[,batch],
                       BPPARAM = parallel_param)
    if(save_plots) {
      p2 <- plotColData(sce, x="scDblFinder.class", y = "sum", colour_by="scDblFinder.class") +
        scale_y_log10() + ggtitle("Doublets")
      ggsave(p2, filename = paste0("./", name, "/doublet_plot.png"), width =  6, height = 6, device = "png")
    }
  }
  return(sce)
}
