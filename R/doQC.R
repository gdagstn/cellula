

#' SCE QC sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param batch character, the name of the column in `colData(sce)` with batch labels
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
                 batch,
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

    if(verbose) cat(blue("[QC/EMPTY]"),"Running emptyDrops. \n")

    if(emptydrops_cutoff == "auto") {
      barcode_ranks <- barcodeRanks(sce)
      emptydrops_cutoff = metadata(barcode_ranks)$inflection
    }

    empty_droplets <- emptyDrops(sce, lower = emptydrops_cutoff)
    keep_droplets <- empty_droplets$FDR <= emptydrops_alpha
    sce$empty <- factor(ifelse(empty_droplets$FDR <= emptydrops_alpha, "ok", "empty"))
    sce$empty[which(is.na(sce$empty))] = "empty"
    if(verbose) print(table(Sig = keep_droplets, Limited = empty_droplets$Limited))
    if(verbose) cat("Empty cells: ", sum(sce$empty == "empty"), "\n")
    sce = sce[, which(sce$empty == "ok"), drop = FALSE]
  }

  if(verbose) cat(blue("[QC]"),"Calculating QC metrics. \n")

  # Generate QC
  if(subset_mito){
    if(sum(grepl("^MT", rowData(sce)$Symbol)) == 0) {
      cat("   No MT genes found.\n")
      mito = NA
    } else mito = rownames(sce)[grepl("^MT-", rowData(sce)$Symbol)]
  } else mito = NA

  if(subset_malat1){
    if(sum(grepl("^MALAT1", rowData(sce)$Symbol)) == 0) {
      cat("   No MALAT1 gene found.\n")
      Malat1 = NA
    } else Malat1 = rownames(sce)[grepl("^MALAT1", rowData(sce)$Symbol)]
  } else Malat1 = NA

  if(subset_ribo){
    if(sum(grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol)) == 0) {
      cat("   No ribo genes found.\n")
      Ribo = NA
    } else Ribo = rownames(sce)[grepl("^MRPL|^MRPS|^RPL|^RPS", rowData(sce)$Symbol)]
  } else Ribo = NA

  subset_list = list(mito = mito, Malat1 = Malat1, Ribo = Ribo)
  subset_list = subset_list[!is.na(subset_list)]

  if(length(subset_list) > 0) {
    sce_fqc <- perCellQCMetrics(sce,
                                subsets = subset_list,
                                BPPARAM = parallel_param)

    low.lib <- isOutlier(log10(sce_fqc$sum),
                         batch = colData(sce)[,batch],
                         type = "lower",
                         nmads=3)

    low.genes <- isOutlier(log10(sce_fqc$detected),
                           batch = colData(sce)[,batch],
                           type = "lower",
                           nmads=3)

    if(!all(is.na(mito))) {
      high.mt <- isOutlier(sce_fqc$subsets_mito_percent,
                           batch = colData(sce)[,batch],
                           type = "higher",
                           nmads = 3)
    } else high.mt <- NA

    if(!all(is.na(Malat1))) {
      high.malat1 <- isOutlier(sce_fqc$subsets_Malat1_percent,
                               batch = colData(sce)[,batch],
                               type = "higher",
                               nmads = 3)
    } else high.malat1 <- NA

    if(!all(is.na(Ribo))) {
      high.ribo <- isOutlier(sce_fqc$subsets_Ribo_percent,
                             batch = colData(sce)[,batch],
                             type = "higher",
                             nmads = 3)
    } else high.ribo <- NA

    data.frame(LowLib=sum(low.lib),
               LowNgenes=sum(low.genes),
               HighMT = ifelse(all(is.na(mito)), FALSE, sum(high.mt)),
               HighMalat1 = ifelse(all(is.na(Malat1)), FALSE, sum(high.malat1)),
               HighRibo = ifelse(all(is.na(Ribo)), FALSE, sum(high.ribo)))

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
  }

  # Doublet finding
  if(detect_doublets){
    if(verbose) cat(blue("[QC/DBL]"), "Finding doublets. \n")

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