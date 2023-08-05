#' SCE Cell Type annotation pipeline
#'
#' Automatic assignment of cell type identities in SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param genesets named list of character vectors, genesets used to calculate
#'     labels and/or assign scores. Alternatively, a character vector to calculate
#'     a single score from.
#' @param method character, the method for scoring. One of \code{"AUC"}, \code{"Seurat"},
#'     \code{"ssGSEA"}, \code{"UCell"}, or \code{"Jaitin"}.
#' @param ref a matrix containing expression values from reference transcriptomes. 
#'     Only used when \code{method = "Jaitin"}.  
#' @param assay the name of the slot in \code{sce} containing expression values 
#'     to calculate the log-likelihood. Only used when \code{method = "Jaitin"} or
#'     \code{method = "ssGSEA"}. 
#' @param verbose logical, should messages on progress be printed? Default is TRUE
#' @param name character, the name of the column in \code{colData(sce)} where final
#'     labels will be stored
#' @param return_scores logical, should the scores for each cell and each geneset
#'     be saved in the object metadata slot? Default is \code{FALSE}.
#' @param kcdf character, which kernel to use for the CDF. One of \code{"Poisson"} or
#'     \code{"Gaussian"}. Only used when \code{method = "ssGSEA"}.
#' @param BPPARAM a \code{BiocParallel} BPPARAM specifying the parallelization. 
#'     Only used when \code{method = "AUC"}. Default is \code{SerialParam()}      
#' @param ... other arguments passed internally to \code{AUCell::AUCell_calcAUC()}
#'     ("AUC" method), \code{Seurat::AddModuleScore()} ("Seurat" method),
#'     \code{GSVA::gsva()}("ssGSEA" method), or \code{UCell::ScoreSignatures_UCell()}
#'     ("UCell" method).
#'
#' @return  a \code{SingleCellExperiment} object with a column named \code{"name"} containing
#'     the highest scoring label for a method. Optionally, the single scores from
#'     each method for each geneset are saved in the \code{metadata} slot.
#'
#' @details This is a wrapper around four methods for assigning to each cell a
#'    score for the expression of a given geneset, and a method that uses a reference
#'    transcriptome as prior probability distribution to estimate the posterior 
#'    probability that a cell is similar to the reference via a multinomial model.
#'    \itemize{
#'     \item{ The "AUC" method
#'     (\href{https://www.nature.com/articles/nmeth.4463}{Aibar et al. 2017}) calculates
#'     the Area Under the Curve of the ranked expression of genes in each geneset.}
#'
#'    \item{The "Seurat" method (\href{https://www.science.org/doi/10.1126/science.aad0501}{Tirosh et al. 2016})
#'     uses the \code{Seurat::AddModuleScore()} function. }
#'
#'    \item{The "ssGSEA" method (\href{https://www.nature.com/articles/nature08460}{Barbie et al. 2009})
#'     calculates an enrichment score in each cell using "single sample" GSEA.}
#'
#'    \item{The "UCell" method (\href{https://www.sciencedirect.com/science/article/pii/S2001037021002816}{Andreatta and Carmona 2020})
#'     calculates an enrichment score using a Mann-Whitney U-test.}
#'     
#'    \item{The "Jaitin" method (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4412462/}{Jaitin et al. 2016}) calculates the posterior probability of a single cell being close to a 
#'     reference transcriptome.}
#'    }
#'
#'    If a list of genesets is given, all methods result in a final label assignment
#'    (i.e. the names of the user-provided geneset for which each single cell has
#'    the highest score). If only a single geneset is provided, the score for that
#'    geneset is added to the \code{colData} slot.
#'    
#'    For the \code{"Jaitin"} method, the \code{ref} object is a matrix containing
#'    genes in each row, and samples in each column. It is important to specify the
#'    \code{assay} slot depending on what the reference contains. If the reference
#'    contains raw counts, it is advisable to use \code{assay = "counts"}; if it 
#'    contains log-normalized counts, it is advisable to use \code{assay = "logcounts"}.
#'    The resulting scores, optionally saved in the \code{metadata}, are the posterior
#'    probabilities. 
#'    
#'    It is possible to use \code{buildReference()} to create a reference matrix
#'    out of a \code{SingleCellExperiment} object to be used as an input.
#'
#'    Please note that scores are not comparable between themselves - AUC,
#'    ssGSEA, UCell and Jaitin are strictly positive, whereas the Module Score can also
#'    be negative; UCell, AUCell and Jaitin are normalized between 0 and 1, whereas
#'    ssGSEA is not.
#'
#'    The user can pass additional arguments to each function using the `...`
#'    argument.
#'
#' @importFrom BiocParallel SerialParam 
#' @export

assignIdentities <- function(sce,
                             genesets = NULL,
                             method,
                             ref = NULL,
                             assay = NULL,
                             verbose = TRUE,
                             name = NULL,
                             return_scores = FALSE,
                             kcdf = "Gaussian",
                             BPPARAM = SerialParam(),
                             ...) {
  # Checks
  ep = .redm("{cellula::assignIdentities()} - ")
  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if (!method %in% c("AUC", "Seurat", "ssGSEA", "UCell", "Jaitin")) 
    stop(ep, "method not recognized - must be one of \"AUC\", \"Seurat\", 
             \"ssGSEA\", \"UCell\", \"Jaitin\"")
  if (is(genesets, "character")) {
    if (is.null(name)) name = "signature"
    genesets = list(genesets)
    names(genesets) = name
  }

  # Assignment module
  switch(method,
         "AUC" = {sce = .assignIdentities_AUC(sce,
                                              genesets,
                                              verbose,
                                              name = name,
                                              return_scores = return_scores,
                                              BPPARAM = BPPARAM,
                                              ...)},
         "Seurat" = {sce = .assignIdentities_Seurat(sce,
                                                    genesets,
                                                    verbose,
                                                    name = name,
                                                    return_scores = return_scores,
                                                    ...)},
         "ssGSEA" = {sce = .assignIdentities_ssGSEA(sce,
                                                    genesets,
                                                    verbose,
                                                    name = name,
                                                    assay = assay,
                                                    return_scores = return_scores,
                                                    kcdf = kcdf,
                                                    ...)},
         "UCell" = {sce = .assignIdentities_UCell(sce,
                                                  genesets,
                                                  verbose,
                                                  name = name,
                                                  return_scores = return_scores,
                                                  ...)},
         "Jaitin" = {sce = .assignIdentities_Jaitin(sce,
                                                   ref,
                                                   assay,
                                                   verbose,
                                                   name = name,
                                                   return_scores = return_scores,
                                                   ...)}
         )

  return(sce)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocParallel SerialParam
#' @importFrom utils tail
#' @importFrom S4Vectors metadata metadata<-

.assignIdentities_AUC <- function(sce,
                                  genesets,
                                  verbose,
                                  name = NULL,
                                  return_scores = FALSE,
                                  BPPARAM = SerialParam(),
                                  ...){
    # Checks
    ep = .redm("{cellula::.assignIdentities_AUC()} - ")
    
    if (!"AUCell" %in% rownames(installed.packages()))
      stop(ep, "the `AUCell` package must be installed first.\n
                Run `BiocManager::install(\"AUCell\") to use this function.")
    if (is.null(name)) labelname = "labels_AUC" else labelname = name
    if (verbose) cat(.bluem("[ANNO/AUC]"), "Assigning cell labels \n")
  
   rankings <- AUCell::AUCell_buildRankings(counts(sce),
                                             splitByBlocks = TRUE,
                                             plotStats = FALSE,
                                             verbose = verbose,
                                             BPPARAM = BPPARAM)
   aucs <- AUCell::AUCell_calcAUC(genesets,
                                  rankings,
                                  aucMaxRank = ceiling(0.2 * nrow(rankings)),
                                  nCores = BPPARAM$workers, 
                                  ...)
   scores <- as.data.frame(t(assay(aucs)))
  
   if (length(genesets) > 1) {
        # Decide scores
        labels_AUC = names(genesets)[apply(scores, 1, which.max)]
        colData(sce)[,labelname] = labels_AUC
       if (return_scores) metadata(sce)$AUC_Scores = scores
      } else if (length(genesets) == 1) {
        colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
      }
      sce
}

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom utils tail
#' @importFrom S4Vectors metadata metadata<-

.assignIdentities_Seurat <- function(sce,
                                     genesets,
                                     verbose,
                                     return_scores = FALSE,
                                     name = NULL,
                                     ...){
  # Checks
  ep = .redm("{cellula::.assignIdentities_Seurat()} - ")
  if(!"Seurat" %in% rownames(installed.packages()))
    stop(ep, "the `Seurat` package must be installed first.\n
                Run `BiocManager::install(\"Seurat\") to use this function.")
  
  if (is.null(name)) labelname = "labels_Seurat" else labelname = name

  if (verbose) cat(.bluem("[ANNO/Seurat]"), "Adding module scores \n")
  old_colnames = colnames(sce)
  if (any(duplicated(colnames(sce)))) {
    colnames(sce) = paste0("cell_", seq_len(ncol(sce)))
  }

  # Conversion to Seurat object
  seu <- SeuratObject::CreateSeuratObject(counts = counts(sce),
                                          meta.data = as.data.frame(colData(sce)))
  logc = logcounts(sce)
  rownames(logc) = rownames(seu)
  seu <- SeuratObject::SetAssayData(object = seu, slot = "data", new.data = logc)
  # Module score
  seu <- Seurat::AddModuleScore(seu, features = genesets, ...)
  colnames(seu@meta.data)[tail(seq_along(seu@meta.data), length(genesets))] = names(genesets)
  scores = seu@meta.data[,tail(seq_along(seu@meta.data), length(genesets))]

  if (length(genesets) > 1) {
    # Decide label
    labels_Seurat = names(genesets)[apply(scores, 1, which.max)]
    colData(sce)[,labelname] = labels_Seurat
    if (return_scores) metadata(sce)$Seurat_ModuleScores = scores
  } else if (length(genesets) == 1) {
    colData(sce)[,labelname] = as.numeric(scores)
  }
  sce
}

#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata metadata<-

.assignIdentities_ssGSEA <- function(sce,
                                     genesets,
                                     verbose,
                                     name = NULL,
                                     assay = "logcounts",
                                     return_scores = FALSE,
                                     kcdf = "Gaussian",
                                     ...) {
  
  # Checks
  ep = .redm("{cellula::.assignIdentities_ssGSEA()} - ")
    if (!"GSVA" %in% rownames(installed.packages()))
    stop(ep, "the `GSVA` package must be installed first.\n
                Run `BiocManager::install(\"GSVA\") to use this function.")
  if (is.null(name)) labelname = "labels_ssGSEA" else labelname = name
  if (is.null(assay))
    stop(ep, "You must specify the assay argument (typically \"logcounts\"")

  if (verbose) cat(.bluem("[ANNO/ssGSEA]"), "Calculating ssGSEA \n")
  ss = GSVA::gsva(sce,
                  gset.idx.list = genesets,
                  annotation = assay,
                  method = "ssgsea",
                  kcdf = kcdf,
                  verbose = verbose,
                  ...)

  scores = t(as.matrix(assay(ss, "es")))

  if (length(genesets) > 1) {
    # Decide labels
    labels_ssGSEA = names(genesets)[apply(scores, 1, which.max)]
    colData(sce)[,labelname] = labels_ssGSEA
    if (return_scores) metadata(sce)$ssGSEA_Scores = scores
  } else if (length(genesets) == 1) {
    colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
  }
  if (verbose) message(.bluem("[ANNO/ssGSEA]"), "Done.")
  sce
}

#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata metadata<-

.assignIdentities_UCell <- function(sce,
                                    genesets,
                                    verbose,
                                    name = NULL,
                                    return_scores = FALSE,
                                    ...) {
  
  ep = .redm("{cellula::.assignIdentities_UCell()} - ")
  
  if (!"UCell" %in% rownames(installed.packages()))
    stop(ep, "the `UCell` package must be installed first.\n
                Run `BiocManager::install(\"UCell\") to use this function.")
  if (is.null(name)) 
    labelname = "labels_UCell" else labelname = name
    
  if (verbose) 
    message(.bluem("[ANNO/UCell]"), "Calculating UCell scores")

  maxrank = 1500
  scores <- UCell::ScoreSignatures_UCell(matrix = assay(sce, "counts"),
                                         features = genesets,
                                         maxRank = maxrank,
                                         name = "",
                                         ...)
  if (length(genesets) > 1) {
    # Decide labels
    labels_UCell = names(genesets)[apply(scores, 1, which.max)]
    colData(sce)[,labelname] = labels_UCell
    if(return_scores) metadata(sce)$UCell_Scores = scores
  } else if (length(genesets) == 1) {
    colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
  }
  if (verbose) message(.bluem("[ANNO/UCell]"), "Done.")
  sce
}

#' @importFrom matrixStats colMaxs
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata

.assignIdentities_Jaitin <- function(sce,
                                     ref,
                                     assay = "counts",
                                     verbose = FALSE,
                                     name = NULL,
                                     return_scores = FALSE){
  #Checks
  ep = .redm("{cellula::assignIdentities()} - ")
  if (is.null(ref)) stop(ep, "a reference must be supplied through the ref argument")
  if (is.null(assay)) stop(ep, "an assay slot must be supplied through the assay argument")
  if(is.null(name)) labelname = "labels_Jaitin" else labelname = name
  
  common = intersect(rownames(ref), rownames(sce))
  if (length(common) == 0) {
    common = intersect(rownames(ref), rowData(sce)$ID)
    if (length(common) == 0) {
      common = intersect(rownames(ref), rowData(sce)$Symbol)
      if (length(common) == 0) {
        stop(ep, "no genes in common were found in any slot (rownames, rowData ID and Symbol)")
      } else {
        old.rownames = rownames(sce)
        rownames(sce) = rowData(sce)$Symbol
      }
    } else {
      old.rownames = rownames(sce)
      rownames(sce) = rowData(sce)$ID
    }
  } else {
    old.rownames = rownames(sce)
  }
  
  if (length(common) < floor(0.1 * nrow(sce))) 
    warning(ep, "less than 10% genes in common between object and reference")
  
  # Scaling the reference
  rl = t(t(ref)/colSums(ref))
  
  # Keep only common genes between sce and reference, set 0 to very low non-neg
  ref_common = rl[common,]
  ref_common[ref_common == 0] = 1e-8
  counts_common = assay(sce, assay)[common,]
  
  if (verbose) cat(.bluem("[ANNO/JAITIN]"), "Calculating log-likelihood \n")

  loglik = as(t(t(counts_common) %*% log(ref_common)), "matrix")
  
  colnames(loglik) = colnames(sce)
  rownames(loglik) = colnames(ref)
  
  if (verbose) message(.bluem("[ANNO/JAITIN]"), "Calculating posterior probabilities")
  
  posterior = exp(t(loglik)-colMaxs(loglik)-log(rowSums(t(loglik)/colMaxs(loglik))))

  if (verbose) message(.bluem("[ANNO/JAITIN]"), "Deciding best labels \n")
  
  toplabel = colnames(posterior)[apply(posterior, 1, which.max)]
  
  if (return_scores) metadata(sce)[[paste0(labelname, "_scores")]] = posterior
  
  if (verbose) message(.bluem("[ANNO/JAITIN]"), "All done \n")
  
  colData(sce)[,labelname] = toplabel
  rownames(sce) = old.rownames
  sce
}

#' Build a reference
#' 
#' Summarizes an assay from a SingleCellExperiment object based on a character column
#' 
#' @param sce a SingleCellExperiment object containing the reference data
#' @param agg_by character, the name of the \code{colData(sce)} column that the 
#'     aggregation should be done by, e.g. a column with labels or clustering
#'     results
#' @param agg_assay character, the name of the \code{assay(sce)} slot whose values
#'     will be aggregated. Default is "logcounts"
#' @param agg_fun character, the aggregation function, one of "sum" or "mean". 
#'     Default is "sum".
#' @returns a matrix with the same row names as \code{sce} and values aggregated by
#'     unique values of \code{agg_by}, which will be the column names.
#'     
#' @importFrom SummarizedExperiment colData assay
#' @importFrom Matrix rowMeans rowSums
#' @export

buildReference <- function(sce, agg_by, agg_assay = "logcounts", agg_fun = "sum"){
    # Checks
    ep = .redm("{cellula::buildReference()} - ")
    if (is.null(agg_by)) 
      stop(ep, "agg_by must be supplied to specify which column to aggregate by")
    if (is.null(agg_assay)) 
      stop(ep, "an assay slot must be supplied through the assay argument")
    if (!agg_by %in% colnames(colData(sce)))
      stop(ep, "agg_by was not found in the colData slot")
    if (!is(colData(sce)[,agg_by], "factor") & !is(colData(sce)[,agg_by], "character"))
      stop(ep, "agg_by must point to either a factor or character in colData")
    if (is(colData(sce)[,agg_by], "factor")) {
      colData(sce)
    }
    # Decide aggregation function
    # calls directly from Matrix to dispatch properly to sparse matrices
    afn = switch(agg_fun,
                 "sum" = Matrix::rowSums,
                 "mean"= Matrix::rowMeans)
    # Loop applying the aggregation function
    quants = SummarizedExperiment::assay(sce, agg_assay)
    labs = unique(as.character(colData(sce)[,agg_by]))
    refmat = do.call(cbind, lapply(labs, function(x) 
           afn(quants[, which(colData(sce)[,agg_by] == x), drop = FALSE])
         ))
  colnames(refmat) = labs
  rownames(refmat) = rownames(sce)
  refmat
}
