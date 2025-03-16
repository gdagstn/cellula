#' SCE QC sub-pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels.
#'     Default is \code{NULL} meaning no batches will be considered, and data will be
#'     processed as a single batch.
#' @param name character, the name of the file/folder.
#' @param discard logical, should values that do not meet QC thresholds be discarded?
#'     Default is \code{TRUE}.
#' @param subset_mito logical, should mitochondrial transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param subset_ribo logical, should ribosomal transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param subset_malat1 logical, should MALAT1 transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param detect_doublets logical, should \code{scDblFinder} be run? Default is \code{TRUE}
#' @param run_emptydrops logical, should `emptyDrops` be run? Default is \code{TRUE}
#' @param emptydrops_cutoff either \code{"auto"} (default, barcode rank inflection point)
#'     or a numeric. Cells with total reads below this cutoff are used to calculate
#'     ambient RNA profiles and are removed.
#' @param emptydrops_alpha numeric, the FDR threshold to call an empty barcode.
#' @param verbose logical, display messages on progress? Default is \code{FALSE}
#' @param save_plots logical, should plots be drawn and saved? Default is \code{TRUE}
#' @param path character, the path to save the plots. Default is the local working directory ("./")
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Default is \code{SerialParam()}
#'     meaning no parallelization will be used.
#'
#' @return  a \code{SingleCellExperiment} object with QC including doublet assignment
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom scuttle isOutlier perCellQCMetrics
#' @importFrom ggplot2 scale_y_log10 ggtitle ggsave ylab
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
				 path = "./",
                 parallel_param = SerialParam()){

  
  ep <- .redm("{cellula::doQC()} - ")
  
  if (detect_doublets){
    dependencies = data.frame("package" = c("scDblFinder"),
                              "repo" = c("BioC"))
    if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))

  }
  
  if (run_emptydrops){
      sce <- doEmptyDrops(sce, batch = batch,
                        emptydrops_cutoff = emptydrops_cutoff,
                        emptydrops_alpha = emptydrops_alpha,
                        verbose = verbose,
                        parallel_param = parallel_param)
  }
  if (verbose) message(.bluem("[QC]")," Calculating QC metrics.")
  # Generate QC
  if (subset_mito){
    if (sum(grepl("^MT-|^Mt-", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      warning("   No MT genes found.\n")
      mito <- FALSE
    } else mito = rownames(sce)[grepl("^MT-|^Mt-", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else mito <- FALSE

  if (subset_malat1){
    if (sum(grepl("^MALAT1|^Malat1", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      warning("   No MALAT1 gene found.\n")
      Malat1 <- FALSE
    } else Malat1 = rownames(sce)[grepl("^MALAT1|^Malat1", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else Malat1 <- FALSE

  if (subset_ribo){
    if (sum(grepl("^MRPL|^MRPS|^RPL|^RPS|^Mrpl|^Mrps|^Rpl|^Rps", rowData(sce)$Symbol, ignore.case = TRUE)) == 0) {
      warning("   No ribo genes found.\n")
      Ribo <- FALSE
    } else Ribo <- rownames(sce)[grepl("^MRPL|^MRPS|^RPL|^RPS|^Mrpl|^Mrps|^Rpl|^Rps", rowData(sce)$Symbol, ignore.case = TRUE)]
  } else Ribo <- FALSE

  subset_list <- list(mito = mito, Malat1 = Malat1, Ribo = Ribo)
  subset_list <- subset_list[!is.na(subset_list)]

  if (length(subset_list) == 0) subset_list = NULL

  if (!is.null(batch)) qcbatch <- colData(sce)[,batch] else qcbatch <- NULL
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
    if (!all(is.na(mito))) {
      high.mt <- isOutlier(sce_fqc$subsets_mito_percent,
                           batch = qcbatch,
                           type = "higher",
                           nmads = 3)
    } else high.mt <- FALSE

    if (!all(is.na(Malat1))) {
      high.malat1 <- isOutlier(sce_fqc$subsets_Malat1_percent,
                               batch = qcbatch,
                               type = "higher",
                               nmads = 3)
    } else high.malat1 <- FALSE

    if (!all(is.na(Ribo))) {
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
    sce$discard[is.na(sce$discard)] = TRUE
    
    # Save plots  
    if (save_plots){
      if (verbose) message(.bluem("[QC]"),"   Saving QC plots.")
      savepath <- paste0(path, "/", name)
      dir.create(paste0(savepath, "/plots"), showWarnings = FALSE, recursive = TRUE)
      savepath <- paste0(savepath, "/plots")
      if (!is.null(batch)) {
        xlen <- length(unique(colData(sce)[,batch])) * 2
        w <- 250 * xlen
        h <- 600 * sqrt(xlen) 
      } else {
        w <- 1200
        h <- 1200
      }
      if (!all(colData(sce)[,"discard"] == FALSE)) cby = "discard" else cby = NULL
      p1 <- plot_Coldata(sce, y = "sum", x = batch, color_by = cby) + 
        scale_y_log10() + ylab("log10(Total UMI)")
      ggsave(filename = "QC_Plots_total_UMI.pdf", p1, 
             path = savepath, device = "pdf",
             width = w, height = h, units = "px")
      p2 <- plot_Coldata(sce, y = "detected", 
                         x = batch, color_by = cby) +
        ylab("Total detected genes")
      ggsave(filename = "QC_Plots_detected_genes.pdf", p2, 
             path = savepath, device = "pdf",
             width = w, height = h, units = "px")
      p3 <- plot_Coldata(sce, y = "sum", x = "detected", group_by = batch) + 
        scale_y_log10() +
        xlab("Total detected genes") + ylab("log10(Total UMI)")
      ggsave(filename = "QC_Plots_detected_by_total_UMI.pdf", p3, 
             path = savepath, device = "pdf",
             width = h*1.1, height = h*1.1, units = "px")
      if (!(all(mito == FALSE))) {
        pmito <- plot_Coldata(sce, y = "subsets_mito_percent", x = batch, 
                              color_by = cby) + 
          ylab("% Mitochondrial transcripts")
        ggsave(filename = "QC_Plots_mito_percent.pdf", pmito, 
               path = savepath, device = "pdf", 
               width = w, height = h, units = "px")
      } else {
        pmito <- NA
      }
      if (!(all(Ribo == FALSE))) {
        pribo <- plot_Coldata(sce, y = "subsets_Ribo_percent", x = batch, 
                              color_by = cby) + 
          ylab("% Ribosomal transcripts")
        ggsave(filename = "QC_Plots_ribo_percent.pdf", pribo, 
               path = savepath, device = "pdf", 
               width = w, height =h, units = "px")
      } else {
        pribo <- NA
      }
      
      if (!(all(Malat1 == FALSE))) {
        pmalat1 <- plot_Coldata(sce, y = "subsets_Malat1_percent", x = batch, 
                                color_by = cby)  + 
          ylab("% MALAT1 transcripts")
        ggsave(filename = "QC_Plots_MALAT1_percent.pdf", pmalat1,
               path = savepath, device = "pdf", 
               width = w, height = h, units = "px")
      } else {
        pmalat1 <- NA
      }
    }

    if(discard) sce <- sce[, !sce$discard]

  # Doublet finding
  if (detect_doublets){
    if (verbose) message(.bluem("[QC/DBL]"), "Finding doublets.")
    if (!is.null(batch)) samples = colData(sce)[,batch] else samples = NULL
    sce <- scDblFinder::scDblFinder(sce, verbose = verbose,
                                    samples = samples,
                                    BPPARAM = parallel_param)
  }
  
     if (save_plots) {
      if (detect_doublets) {
        pdbl <- plot_Coldata(sce, y = "sum", x = batch, color_by = "scDblFinder.class") + 
          scale_y_log10()  + 
          ylab("log10(Total UMI)")
        ggsave(filename = "QC_Plots_total_UMI_doublets.pdf", pdbl, 
               path = savepath, device = "pdf", 
               width = w, height = h, units = "px")
      } else {
        pdbl <- NA
      }
      pl <- list(p1, p2, p3, pmito, pribo, pmalat1, pdbl)
      pl <- pl[!is.na(pl)]
      ncols <- floor(sqrt(length(pl)))
      pfinal <- do.call(arrangeGrob, c(pl, ncol = ncols))
      ggsave(filename = "QC_Plot_ALL.pdf", pfinal, 
             path = savepath, device = "pdf", 
             width = 3200, height = 3200, units = "px")
  } 
  
  sce
}
