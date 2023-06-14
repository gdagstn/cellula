#'#' SCE Cell Type annotation pipeline
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
#'     to calculate the log-likelihood. Only used when \code{method = "Jaitin"}.
#' @param verbose logical, should messages on progress be printed? Default is TRUE
#' @param name character, the name of the column in \code{colData(sce)} where final
#'     labels will be stored
#' @param return_scores logical, should the scores for each cell and each geneset
#'     be saved in the object metadata slot? Default is \code{FALSE}.
#' @param kcdf character, which kernel to use for the CDF. One of \code{"Poisson"} or
#'     \code{"Gaussian"}. Only used when \code{method = "ssGSEA"}.
#' @param annotation character, which assay name to use for rank calculation.
#'     Only used when \code{method = "ssGSEA"}.
#' @param BPPARAM a \code{BPPARAM} object from \code{BiocParallel}. Default is 
#'     \code{SerialParam()}, meaning no parallelization will be used.
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
#'    Please note that scores are not comparable between themselves - AUC,
#'    ssGSEA, UCell and Jaitin are strictly positive, whereas the Module Score can also
#'    be negative; UCell, AUCell and Jaitin are normalized between 0 and 1, whereas
#'    ssGSEA is not.
#'
#'    The user can pass additional arguments to each function using the `...`
#'    argument.
#'
#' @importFrom BiocParallel SerialParam BPPARAM
#' @export

assignIdentities <- function(sce,
                             genesets,
                             method,
                             ref = NULL,
                             assay = NULL,
                             verbose = TRUE,
                             name = NULL,
                             return_scores = FALSE,
                             kcdf = "Gaussian",
                             annotation = "logcounts",
                             BPPARAM = SerialParam(),
                             ...) {
  ## Sanity checks
  # Error prefix
  ep = "{cellula::assignIdentities()} - "
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(method %in% c("AUC", "Seurat", "ssGSEA", "UCell", "Jaitin")) 
    stop(paste0(ep, "method not recognized - must be one of \"AUC\", \"Seurat\", \"ssGSEA\", \"UCell\", \"Jaitin\""))

  
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
                                                  ...)},
         "Jaitin" = {sce = .assignIdentities_Jaitin(sce,
                                                   ref,
                                                   assay,
                                                   verbose,
                                                   name = name,
                                                   return_scores = return_scores,
                                                   BPPARAM = BPPARAM
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

#' @importFrom BiocParallel SerialParam MulticoreParam bplapply
#' @importFrom crayon blue
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors metadata

.assignIdentities_Jaitin <- function(sce,
                                     ref,
                                     assay = "counts",
                                     verbose = FALSE,
                                     name = NULL,
                                     return_scores = FALSE,
                                     BPPARAM = SerialParam()){
  
  if(is.null(name)) labelname = "labels_Jaitin" else labelname = name
  
  ep = "{cellula::assignIdentities()} - "
  
  if(is.null(ref)) stop(paste0(ep, "a reference must be supplied through the ref argument"))
  if(is.null(assay)) stop(paste0(ep, "an assay slot must be supplied through the assay argument"))
  
  common = intersect(rownames(ref), rownames(sce))
  
  if(length(common) == 0) {
    common = intersect(rownames(ref), rowData(sce)$ID)
    if(length(common) == 0) {
      common = intersect(rownames(ref), rowData(sce)$Symbol)
      if(length(common) == 0) stop(paste0(ep, "no genes in common were found in any slot (rownames, rowData ID and Symbol)"))
    }
  }
  
  if(length(common) < floor(0.1 * nrow(sce))) warning(paste0(ep, "less than 10% genes in common between object and reference"))
  
  ref_common = ref[common,]
  counts_common = assay(sce, assay)[common,]
  
  if(verbose) cat(blue("[ANNO/JAITIN]"), "Calculating log-likelihood \n")
  
  rl = log(t(t(ref_common) / colSums(ref_common))+1)
  
  loglik = bplapply(seq_len(ncol(sce)), function(x) 
    colSums(counts_common[,x] * rl), 
    BPPARAM = BPPARAM)
  
  if(verbose) cat(blue("[ANNO/JAITIN]"), "Calculating posterior probabilities \n")
  
  posterior = bplapply(loglik, function(x) 
    exp(x-max(x)-log(sum(x/max(x)))), 
    BPPARAM = BPPARAM)
  
  cmat = do.call(rbind, posterior)
  
  if(verbose) cat(blue("[ANNO/JAITIN]"), "Deciding best labels \n")
  
  cmat_top = colnames(cmat)[apply(cmat, 1, which.max)]
  
  if(return_scores) metadata(sce)$JaitinScores = cmat
  
  if(verbose) cat(blue("[ANNO/JAITIN]"), "All done \n")
  
  colData(sce)[,labelname] = cmat_top
  
  return(sce)
}