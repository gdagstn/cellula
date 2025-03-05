#' SingleCellExperiment processing pipeline
#'
#' Pipeline for automatic processing and integration of SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param batch character, the name of the column in \code{colData(sce)} with batch labels
#' @param name character, the name of the project/folder where files will be saved.
#'     If \code{NULL}, a random name will be generated
#' @param do_qc logical, should QC steps be performed? Default is \code{TRUE}
#' @param do_norm logical, should normalization and dimensionality reduction be
#'     performed? Default is \code{TRUE}
#' @param geneset_list named list of gene IDs (must coincide with rownames of
#'     \code{sce}) that will be used by \code{AUCell}. Default is \code{NULL} meaning no AUC
#'     will be calculated.
#' @param discard logical, should values that do not meet QC thresholds be discarded?
#'     Default is \code{TRUE}.
#' @param subset_mito logical, should mitochondrial transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param subset_ribo logical, should ribosomal transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param subset_malat1 logical, should MALAT1 transcripts be used for QC?
#'     Default is \code{TRUE}
#' @param detect_doublets logical, should \code{scDblFinder} be run? Default is \code{TRUE}
#' @param run_emptydrops logical, should \code{emptyDrops} be run? Default is \code{FALSE}.
#' @param emptydrops_cutoff either \code{"auto"} (default, barcode rank inflection point)
#'     or a numeric. Cells with total reads below this cutoff are used to calculate
#'     ambient RNA profiles and are removed.
#' @param emptydrops_alpha numeric, the FDR threshold to call an empty barcode.
#' @param hvg_ntop numeric, the number of top highly variable genes to be selected.
#'     Default is 2000.
#' @param integration_method character, one of \code{"fastMNN"}, \code{"Harmony"}, \code{"Seurat"},
#'     \code{"LIGER"}, or \code{"regression"}. See \code{?integrateSCE} for more 
#'     information.
#' @param ndims numeric, the number of dimensions to retain in the reduced dimension
#'     embedding for downstream applications. Default is 20.
#' @param verbose logical, display messages on progress? Default is \code{FALSE}.
#' @param save_plots logical, should plots be drawn and saved? Default is \code{TRUE}
#' @param save_temp logical, should temporary files be saved as the pipeline progresses?
#'     Default is \code{FALSE}
#' @param path character, the path where the files will be saved. Default is the current
#'     working directory.
#' @param parallel_param a \code{BiocParallel} object specifying the parallelization backend
#'     to be used in some steps of the pipeline. Default is \code{SerialParam()}
#'     meaning no parallelization will be used. Note: for Seurat options, the
#'     \code{future} framework should be set up with maximum size and number of cores.
#'
#' @return  a \code{SingleCellExperiment} object with normalized data, doublet assignment
#'    (if calculated), uncorrected and corrected PCA and UMAP coordinates according
#'    to the method of choice.
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom scuttle isOutlier perCellQCMetrics
#' @importFrom gridExtra grid.arrange
#' @importFrom scater runPCA 
#' @importFrom uwot umap
#' @importFrom ggplot2 scale_y_log10 ggtitle ggsave
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom BiocParallel SerialParam
#'
#' @export

cellula <- function(sce,
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
                     save_temp = FALSE,
					 path = "./",
                     parallel_param = SerialParam()) {

  # Checks (mostly delegated to single modules)
    if (!is.null(batch)) {
   if(!batch %in% colnames(colData(sce)))
    stop("Batch label \"", batch, "\" not found.")
   if (!is.factor(colData(sce)[,batch]))
     colData(sce)[,batch] <- as.factor(colData(sce)[,batch])
 }
  
  # Start parameter logging - not fully implemented
  if(is.null(metadata(sce)$cellula_log)) {
    clog <- .initLog()
  } else clog <- metadata(sce)$cellula_log
  
  # Make folder
  if(is.null(name)) {
    name <- .randomName()
    message("No name selected, so the randomly assigned name is: ", name)
  }

  dir.create(name)
  clog$name <- name
  clog$dir <- paste0(getwd(), "/name")
  
  # Begin
    if(verbose) {
    message("\nWorking on object ", name, "\n")
    ncells <- ncol(sce)
    message("Input cells: ", ncells, "\n")
      clog$qc$input_cells <- ncells
    if(!is.null(batch)) {
      message("By batch: \n")
      print(table(colData(sce)[,batch]))
      clog$qc$input_by_batch <- table(colData(sce)[,batch])
    }
  }
 
  # QC module
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
				path = path,
                parallel_param = parallel_param)
    
    clog$qc$do_qc <- do_qc
    clog$qc$discard <- discard
    clog$qc$subset_mito <- subset_mito
    clog$qc$subset_ribo <- subset_ribo
    clog$qc$subset_malat1 <- subset_malat1
    clog$qc$detect_doublets <- detect_doublets
    clog$qc$run_emptydrops <- run_emptydrops
    clog$qc$emptydrops_cutoff <- emptydrops_cutoff
    clog$qc$emptydrops_alpha <- emptydrops_alpha
    clog$qc$parallel_param <- parallel_param
    clog$qc$output_cells <- ncol(sce)
    if(!is.null(batch)) clog$qc$output_by_batch <- table(colData(sce)[,batch])
    
    metadata(sce)$cellula_log <- clog
    
    if(save_temp) {
      if(verbose) message("Saving temporary file.")
      saveRDS(sce, file = paste0(path, "/", name, "/", name, "_tempSCE.RDS"))
    }
  }

  # Normalization and dimensionality reduction module
    if(do_norm) {
    sce <- doNormAndReduce(sce, batch = batch, name = name,
                           ndims = ndims,
                           hvg_ntop = hvg_ntop,
                           verbose = verbose,
                           parallel_param = parallel_param)

    if(save_temp) {
      if(verbose) message("Saving temporary file. \n")
      saveRDS(sce, file = paste0(path, "/", name, "/", name, "_tempSCE.RDS"))
    }
  }
  hvgs <- metadata(sce)$hvgs
  
  # Integration module
    if(!is.null(batch)) {
    if(verbose) message(.bluem("[INT] "), "Integration.")
    sce <- integrateSCE(sce,
                       batch = batch,
                       hvgs = hvgs,
                       hvg_ntop = hvg_ntop,
                       method = integration_method,
                       ndims = ndims,
                       parallel_param = parallel_param,
                       verbose = verbose)
  }

  if(verbose) message("Saving final object.")
  	saveRDS(sce, file = paste0(path, "/", name, "/", name, "_PS_INT_SCE.RDS"))

  if (save_temp) {
    if (verbose) message("Deleting temporary file.")
    file.remove(paste0(path, "/", name, "/", name, "_tempSCE.RDS"))
  }

  if (verbose) message("All done. Input cells: ", ncells, ", final cell number: ", ncol(sce))
  metadata(sce)$cellula_log[["output_cells"]] <- ncol(sce)
  if (!is.null(batch)) {
    metadata(sce)$cellula_log[["output_by_batch"]] <- table(colData(sce)[,batch])
  }
  sce
}

#' @noRd
#' @importFrom utils sessionInfo
.initLog <- function() {
  
  list("name" = NULL,
       "dir" = NULL,
       "batch" = NULL,
       
       "qc" <- list("do_qc" = NULL,
                   "input_cells" = NULL,
                   "input_by_batch" = NULL,
                   "output_cells" = NULL,
                   "output_by_batch" = NULL,
                   "discard" = NULL,
                   "subset_mito" = NULL,
                   "subset_ribo" = NULL,
                   "subset_malat1" = NULL,
                   "detect_doublets" = NULL,
                   "run_emptydrops" = NULL,
                   "emptydrops_cutoff" = NULL,
                   "emptydrops_alpha" = NULL,
                   "save_plots" = NULL),
    
       
       "norm_reduce" <- list("hvg_ntop" = NULL,
                            "pca" = NULL,
                            "umap_min_dist" = NULL,
                            "umap_n_neighbors" = NULL,
                            "umap_other" = NULL,
                            "parallel_param" = NULL),
       
      
       "integration" <- list("integration_method" = NULL,
                            "ndims" = NULL,
                            "parallel_param" = NULL),
       
       "clustering" <- list( "neighbors" = NULL,
                            "weighting_scheme" = NULL,
                            "clustering_sweep" = NULL,
                            "graph_ks" = NULL,
                            "resolution_ks" = NULL,
                            "clustering_method" = NULL,
                            "calculate_modularity" = NULL,
                            "calculate_silhouette" = NULL,
                            "leiden_iterations" = NULL,
                            "save_graphs" = NULL,
                            "prefix" = NULL,
                            "metacluster" = NULL,
                            "metacluster_clusters" = NULL,
                            "metacluster_threshold" = NULL,
                            "metacluster_denominator" = NULL,
                            "parallel_param" = NULL),
       
       "assign_identities" <- list("method" = NULL,
                                  "genesets" = NULL,
                                  "ref" = NULL,
                                  "assay" = NULL,
                                  "name" = NULL,
                                  "return_scores" = NULL,
                                  "kcdf" = NULL,
                                  "parallel_param" = NULL),
       
       "trajectories" <- list("method" = NULL,
                             "dr" = NULL,
                             "clusters" = NULL,
                             "ndims" = NULL,
                             "dr_embed" = NULL,
                             "start" = NULL,
                             "Monocle_lg_control" = NULL,
                             "omega" = NULL,
                             "omega_scale" = NULL,
                             "do_de" = NULL,
                             "de_method" = NULL,
                             "batch_de" = NULL,
                             "parallel_param" = NULL),
  
       "sessionInfo" <- sessionInfo()
  )
}
