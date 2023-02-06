#'#' SCE Cell Type annotation pipeline
#'
#' Automatic assignment of cell type identities in SingleCellExperiment objects
#'
#' @param sce a SingleCellExperiment object
#' @param genesets named list of character vectors, genesets used to calculate 
#'     labels and/or assign scores
#' @param method character, the method for label calling. One of "AUC" or "Seurat".
#' @param verbose logical, should messages on progress be printed? Default is TRUE
#' @param name character, the name of the column in the colData of sce where final
#'     labels will be stored
#' @param ... other arguments passed internally to \code{AUCell::AUCell_calcAUC()} 
#'     ("AUC" method) or \code{Seurat::AddModuleScore()} ("Seurat" method)
#'     
#' @return  a `SingleCellExperiment` object with a column named `"name"` containing
#'     the highest scoring label for a method, and other method-specific columns
#'     such as AUC values or score values per geneset.
#'     
#' @details This is a wrapper around two method for assigning to each cell a 
#'    score for the expression of a given geneset. The "AUC" method calculates
#'    the Area Under the Curve of the ranked expression of genes in each geneset,
#'    whereas the "Seurat" function uses the \code{Seurat::AddModuleScore()}
#'    function. Both methods result in a final label assignment (i.e. the names
#'    of the user-provided geneset for which each single cell has the highest 
#'    score), although scores are not comparable between themselves - AUC is
#'    strictly positive, whereas the Module Score can also be negative. 
#'    The user can pass additional arguments to each function using the `...` 
#'    argument.
#'     
#' @export     

assignIdentities <- function(sce, 
                             genesets, 
                             method, 
                             verbose = TRUE,
                             name = NULL,
                             ...) {
  
  switch(method, 
         "AUC" = {sce = .assignIdentities_AUC(sce, 
                                              genesets, 
                                              verbose, 
                                              name,
                                              ...)},
         "Seurat" = {sce = .assignIdentities_Seurat(sce, 
                                                    genesets, 
                                                    verbose, 
                                                    name,
                                                    ...)})
  return(sce)
}
  
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom SummarizedExperiment colData
#' @importFrom utils tail

.assignIdentities_AUC <- function(sce, 
                                  genesets, 
                                  verbose, 
                                  name = NULL, 
                                  AUC_threshold = 0.2, ...){
  
    if(is.null(name)) labelname = "labels_AUC"
    if(verbose) cat(blue("[ANNO/AUC]"), "Assigning cell labels \n")
    
    rankings <- AUCell_buildRankings(counts(sce),
                                     plotStats = FALSE,
                                     verbose = FALSE)
    
    aucs <- AUCell_calcAUC(genesets,
                           rankings,
                           aucMaxRank = ceiling(0.2 * nrow(rankings)),
                           ...)
    # All assignments
    
    assigned <- as.data.frame(t(assay(aucs)))
    
    
    # Best overall score
    assigned$first_max_score <- apply(assigned, 1, max)
    
    # Second best score
    assigned$second_max_score <- apply(assigned[,1:(ncol(assigned) - 1)], 1, function(x) {
      return(max(x[x != max(x)]))
    })
    
    # Ambiguous labels
    assigned$ambiguous <- (assigned$first_max_score - assigned$second_max_score)/(assigned$first_max_score + assigned$second_max_score) <= AUC_threshold
    
    assigned$best_label = colnames(assigned)[apply(assigned[,1:(ncol(assigned) - 3)], 1, which.max)]
    
    colData(sce)[,labelname] <- factor(assigned$best_label)
    
    return(sce)
}

#' @importFrom SeuratObject CreateSeuratObject SetAssayData
#' @importFrom Seurat AddModuleScore 
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom utils tail

.assignIdentities_Seurat <- function(sce, genesets, verbose, name = NULL, ...){
  
  if(is.null(name)) labelname = "labels_Seurat"
  
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
  
  colData(sce) = cbind(colData(sce), seu@meta.data[,tail(seq_along(seu@meta.data), length(genesets))])
  colData(sce)[,labelname] = names(genesets)[apply(colData(sce)[,tail(seq_along(colData(sce)),  length(genesets))],
                             1, which.max)]
  
  
  return(sce)
}