#'#' SCE Cell Type annotation pipeline
#'
#' Automatic assignment of cell type identities in SingleCellExperiment objects
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param genesets named list of character vectors, genesets used to calculate
#'     labels and/or assign scores. Alternatively, a character vector to calculate
#'     a single score from.
#' @param method character, the method for scoring. One of \code{"AUC"}, \code{"Seurat"},
#'     \code{"ssGSEA"} or \code{"UCell"}.
#' @param verbose logical, should messages on progress be printed? Default is TRUE
#' @param name character, the name of the column in \code{colData(sce)} where final
#'     labels will be stored
#' @param return_scores logical, should the scores for each cell and each geneset
#'     be saved in the object metadata slot? Default is \code{FALSE}.
#' @param kcdf character, which kernel to use for the CDF. One of \code{"Poisson"} or
#'     \code{"Gaussian"}. Only used when \code{method = "ssGSEA"}.
#' @param annotation character, which assay name to use for rank calculation.
#'     Only used when \code{method = "ssGSEA"}.
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
#'    score for the expression of a given geneset.
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
#'    }
#'
#'    If a list of genesets is given, all methods result in a final label assignment
#'    (i.e. the names of the user-provided geneset for which each single cell has
#'    the highest score). If only a single geneset is provided, the score for that
#'    geneset is added to the \code{colData} slot.
#'
#'    Please note that scores are not comparable between themselves - AUC,
#'    ssGSEA, and UCell are strictly positive, whereas the Module Score can also
#'    be negative; UCell and AUCell are normalized between 0 and 1, whereas
#'    ssGSEA is not.
#'
#'    The user can pass additional arguments to each function using the `...`
#'    argument.
#'
#'
#' @export

assignIdentities <- function(sce,
                             genesets,
                             method,
                             verbose = TRUE,
                             name = NULL,
                             return_scores = FALSE,
                             kcdf = "Gaussian",
                             annotation = "logcounts",
                             ...) {
  ## Sanity checks
  # Error prefix
  ep = "{cellula::assignIdentities()} - "
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(method %in% c("AUC", "Seurat", "ssGSEA", "UCell")) 
    stop(paste0(ep, "method not recognized - must be one of \"AUC\", \"Seurat\", \"ssGSEA\", \"UCell\""))

  
  if(is(genesets, "character")) {
    if(is.null(name)) name = "signature"
    genesets = list(genesets)
    names(genesets) = name
  }

  switch(method,
         "AUC" = {sce = .assignIdentities_AUC(sce,
                                              genesets,
                                              verbose,
                                              name = name,
                                              return_scores = return_scores,
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
                                                    return_scores = return_scores,
                                                    annotation = annotation,
                                                    kcdf = kcdf,
                                                    ...)},
         "UCell" = {sce = .assignIdentities_UCell(sce,
                                                  genesets,
                                                  verbose,
                                                  name = name,
                                                  return_scores = return_scores,
                                                  ...)}
         )

  return(sce)
}

#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail
#' @importFrom S4Vectors metadata

.assignIdentities_AUC <- function(sce,
                                  genesets,
                                  verbose,
                                  name = NULL,
                                  return_scores = FALSE,
                                  ...){

    if(is.null(name)) labelname = "labels_AUC" else labelname = name
    if(verbose) cat(blue("[ANNO/AUC]"), "Assigning cell labels \n")

    rankings <- AUCell_buildRankings(counts(sce),
                                     splitByBlocks = TRUE,
                                     plotStats = FALSE,
                                     verbose = verbose)

    aucs <- AUCell_calcAUC(genesets,
                           rankings,
                           aucMaxRank = ceiling(0.2 * nrow(rankings)),
                           ...)

    scores <- as.data.frame(t(assay(aucs)))

    if(length(genesets) > 1) {

      labels_AUC = names(genesets)[apply(scores, 1, which.max)]

      colData(sce)[,labelname] = labels_AUC

      if(return_scores) metadata(sce)$AUC_Scores = scores

    } else if(length(genesets) == 1) {

      colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
    }

    return(sce)
}

#' @importFrom SeuratObject CreateSeuratObject SetAssayData
#' @importFrom Seurat AddModuleScore
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom utils tail
#' @importFrom crayon blue
#' @importFrom S4Vectors metadata

.assignIdentities_Seurat <- function(sce,
                                     genesets,
                                     verbose,
                                     return_scores = FALSE,
                                     name = NULL,
                                     ...){

  if(is.null(name)) labelname = "labels_Seurat" else labelname = name

  if(verbose) cat(blue("[ANNO/Seurat]"), "Adding module scores \n")

  old_colnames = colnames(sce)
  if(any(duplicated(colnames(sce)))) {
    colnames(sce) = paste0("cell_", seq_len(ncol(sce)))
  }

  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  logc = logcounts(sce)
  rownames(logc) = rownames(seu)
  seu <- SetAssayData(object = seu, slot = "data", new.data = logc)
  seu <- AddModuleScore(seu, features = genesets, ...)

  colnames(seu@meta.data)[tail(seq_along(seu@meta.data), length(genesets))] = names(genesets)

  scores = seu@meta.data[,tail(seq_along(seu@meta.data), length(genesets))]

  if(length(genesets) > 1) {

    labels_Seurat = names(genesets)[apply(scores, 1, which.max)]

    colData(sce)[,labelname] = labels_Seurat

    if(return_scores) metadata(sce)$Seurat_ModuleScores = scores

  } else if(length(genesets) == 1) {

    colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
  }

  return(sce)
}

#' @importFrom GSVA gsva
#' @importFrom crayon blue
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata

.assignIdentities_ssGSEA <- function(sce,
                                     genesets,
                                     verbose,
                                     name = NULL,
                                     annotation = "logcounts",
                                     return_scores = FALSE,
                                     kcdf = "Gaussian",
                                     ...) {

  if(is.null(name)) labelname = "labels_ssGSEA" else labelname = name

  if(verbose) cat(blue("[ANNO/ssGSEA]"), "Calculating ssGSEA \n")

  ss = gsva(sce,
            gset.idx.list = genesets,
            annotation = annotation,
            method = "ssgsea",
            kcdf = kcdf,
            verbose = verbose,
            ...)

  scores = t(as.matrix(assay(ss, "es")))

  if(length(genesets) > 1) {

    labels_ssGSEA = names(genesets)[apply(scores, 1, which.max)]

    colData(sce)[,labelname] = labels_ssGSEA

    if(return_scores) metadata(sce)$ssGSEA_Scores = scores

  } else if(length(genesets) == 1) {

    colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
  }

  if(verbose) cat(blue("[ANNO/ssGSEA]"), "Done. \n")

  return(sce)

}

#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom crayon blue
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata

.assignIdentities_UCell <- function(sce,
                                    genesets,
                                    verbose,
                                    name = NULL,
                                    return_scores = FALSE,
                                    ...) {

  if(is.null(name)) labelname = "labels_UCell" else labelname = name

  if(verbose) cat(blue("[ANNO/UCell]"), "Calculating UCell scores \n")

  maxrank = 1500

  scores <- ScoreSignatures_UCell(matrix = assay(sce, "counts"),
                                  features = genesets,
                                  maxRank = maxrank,
                                  name = "",
                                  ...)

  if(length(genesets) > 1) {

    labels_UCell = names(genesets)[apply(scores, 1, which.max)]

    colData(sce)[,labelname] = labels_UCell

    if(return_scores) metadata(sce)$UCell_Scores = scores

  } else if(length(genesets) == 1) {

    colData(sce)[,labelname] = as.numeric(scores[, 1, drop = TRUE])
  }

  if(verbose) cat(blue("[ANNO/UCell]"), "Done. \n")

  return(sce)
}
