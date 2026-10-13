#' Pseudo-bulk differential expression analysis
#'
#' @description This function performs differential expression
#'     analysis on a single-cell experiment with different replicates,
#'     aggregated into pseudobulk on a per-cluster basis.
#'
#' @param sce a SingleCellExperiment object
#' @param replicates character, column in \code{colData(sce)}
#'     that contains the replicate information
#' @param labels character, column in \code{colData(sce)}
#'     that contains the label information
#' @param condition character, column in \code{colData(sce)}
#'     that contains the condition information
#' @param min_cells numeric, minimum number of cells in a pseudobulk sample
#'     (one replicate in one condition and label) for the sample to be kept.
#'     Samples with fewer cells are removed before testing. Default is 10.
#' @param min_rep numeric, minimum number of replicates per condition and label.
#'     Labels with fewer replicates in either of the conditions of the contrast
#'     are removed. Default is 3.
#' @param subset_conditions logical, should the object be subset to contain
#'     only the conditions specified in the contrast?
#'     Default is \code{FALSE}.
#' @param contrast character, the contrast to be tested.
#'     Default is the second level of the condition minus the first level.
#'     See Details for how it is formatted.
#' @param design a design formula for the differential expression analysis.
#'     Default is NULL, meaning the design will be created automatically
#'     with an intercept term set to 0, e.g. \code{~0 + condition}.
#' @param verbose logical, display messages on progress?
#'     Default is \code{TRUE}.
#' @param parallel_param a \code{BiocParallel} object specifying
#'     the parallelization backend
#'
#' @return a list of data frames with the differential expression results
#'     for each label, obtained through the \code{\link[scran]{pseudoBulkDGE}}
#'     function from \code{{scran}}. The \code{conditions} and \code{inputs}
#'     attributes store the conditions used in the contrast and the names of
#'     the replicates, labels and condition columns.
#'
#' @details This function is a wrapper around the \code{pseudoBulkDGE()}
#'    function from \code{{scran}}. It requires the user to specify the
#'    columns in the \code{colData(sce)} that contain the replicate
#'    information, the label information, and the condition information.
#'
#'    The aggregation happens by taking the combination of the condition,
#'    replicate, and label columns as individual pseudobulk samples. The
#'    differential expression analysis is then performed on these samples
#'    assuming \code{condition} is the categorical variable to use for DE.
#'    Any number of covariates can be included in the \code{design} argument.
#'
#'    Pseudobulk samples with fewer than \code{min_cells} cells are removed.
#'    Labels with fewer than \code{min_rep} samples in either condition of the
#'    contrast are then removed, and so are the samples of other conditions
#'    whose group (label and condition) has fewer than \code{min_rep} samples.
#'    Unless \code{subset_conditions = TRUE}, a condition that is absent from
#'    at least one label is removed from all labels, unless it is part of the
#'    contrast.
#'
#'    An option is given to subset the object to only contain the two conditions
#'    included in the contrast. The difference between subsetting and not
#'    lies in the fact that dispersion estimates will only be calculated for
#'    the samples included in the subset object, and not in the whole object.
#'    This can result in different results (e.g. fewer or more DE genes).
#'    According to current best practices it is NOT recommended to subset
#'    the object, but the option is given.
#'
#'    The contrast to be tested can be specified as well, in the format
#'    \code{"conditionLevel1-conditionLevel2"} where condition is the name of
#'    the column containing the condition, and Level1 and Level2 are the two
#'    levels to compare (levels containing \code{"-"} are not supported). A
#'    single term (\code{"conditionLevel1"}) is taken to be compared against the
#'    first level of the condition. If the contrast is not provided, the second
#'    level of the condition is tested against the first: as
#'    \code{"conditionLevel2-conditionLevel1"} for designs without an intercept,
#'    and as the coefficient \code{"conditionLevel2"} for designs with one.
#'
#'    The default design is a formula with an intercept term set to 0, e.g.
#'    \code{~0 + condition} as it forces the user to explicitly declare in the
#'    contrast which is the reference level. However, if design is provided,
#'    the user can specify a design with an intercept (e.g. \code{~condition})
#'    which means that the reference level for the factor will be taken as the
#'    reference level for the contrast too. Not knowing what the reference
#'    value for the factor is can cause the function to fail if the reference
#'    level is used in the \code{contrast} argument.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment sizeFactors
#' @importFrom scran pseudoBulkDGE
#' @importFrom scuttle computeLibraryFactors logNormCounts aggregateAcrossCells
#' @importFrom scuttle librarySizeFactors
#' @importFrom S4Vectors metadata
#' @importFrom BiocParallel SerialParam
#' @importFrom stats formula terms model.matrix
#' @importFrom methods is
#'
#' @export

doPBDGE <- function(sce,
                    replicates,
                    labels,
                    condition,
                    min_cells = 10,
                    min_rep = 3,
                    subset_conditions = FALSE,
                    contrast = NULL,
                    design = NULL,
                    verbose = TRUE,
                    parallel_param = SerialParam()) {
    ## Sanity checks
    # Error prefix
    ep <- .redm("cellula::doPBDGE() - ")

    # Checks
    if (!is(sce, "SingleCellExperiment")) {
        stop(paste0(ep, "Must provide a SingleCellExperiment object."))
    }
    if (!is.character(replicates)) {
        stop(paste0(ep, "replicates must be a character."))
    }
    if (!is.character(labels)) {
        stop(paste0(ep, "labels must be a character."))
    }
    if (!is.character(condition)) {
        stop(paste0(ep, "condition must be a character."))
    }
    if (!condition %in% colnames(colData(sce))) {
        stop(paste0(ep, "The condition column was not found in the colData."))
    }
    if (!replicates %in% colnames(colData(sce))) {
        stop(paste0(ep, "The replicates column was not found in the colData."))
    }
    if (!labels %in% colnames(colData(sce))) {
        stop(paste0(ep, "The labels column was not found in the colData."))
    }
    if (!is.null(contrast) && !is.character(contrast)) {
        stop(paste0(ep, "The contrast must be a character."))
    }
    if (!is.null(design) && !inherits(design, "formula")) {
        stop(paste0(ep, "The design must be a formula."))
    }
    if (!is.numeric(min_cells) || !is.numeric(min_rep)) {
        stop(paste0(ep, "min_cells and min_rep must be numeric."))
    }

    # Condition levels, design and contrast
    cond_levels <- levels(factor(colData(sce)[, condition]))
    if (length(cond_levels) < 2) {
        stop(paste0(ep, "The condition column must contain at least two levels."))
    }

    if (is.null(design)) {
        design <- formula(paste0("~0 + ", condition))
    }

    if (is.null(contrast)) {
        if (attr(terms(design), "intercept") == 1) {
            contrast <- paste0(condition, cond_levels[2])
        } else {
            contrast <- paste0(condition, cond_levels[2], "-", condition, cond_levels[1])
        }
        if (verbose) message(.bluem("[DE] "), "No contrast provided, testing: ", contrast)
    }

    contrast_terms <- trimws(strsplit(contrast, "-", fixed = TRUE)[[1]])
    conditions_used <- substring(contrast_terms, nchar(condition) + 1)
    if (length(conditions_used) < 1 || length(conditions_used) > 2 ||
        !all(startsWith(contrast_terms, condition)) || !all(conditions_used %in% cond_levels)) {
        stop(paste0(ep, "The contrast must be formatted as \"", condition, "Level1-", condition,
                    "Level2\" (or \"", condition, "Level1\" to compare with the first level), where the levels ",
                    "are levels of the condition column and do not contain \"-\"."))
    }
    if (length(conditions_used) == 1) conditions_used <- c(conditions_used, cond_levels[1])

    if (subset_conditions) {
        if (verbose) {
            message(.bluem("[DE] "),
                    "The object will be subset to contain the following conditions only: \n",
                    paste0("\t", paste(conditions_used, collapse = ", ")))
        }
        sce <- sce[, colData(sce)[, condition] %in% conditions_used]
        cond_levels <- intersect(cond_levels, conditions_used)
    }

    ## Aggregation
    if (verbose) message(.bluem("[DE] "), "Aggregating into pseudobulk.")

    agg <- aggregateAcrossCells(sce,
        ids = paste0(colData(sce)[, condition], "__",
                     colData(sce)[, replicates], "__",
                     colData(sce)[, labels]),
        statistics = "sum",
        use.assay.type = "counts",
        BPPARAM = parallel_param
    )

    colData(agg)[[condition]] <- factor(agg[[condition]], levels = cond_levels)
    colData(agg)[[labels]] <- factor(agg[[labels]])

    ## Filtering
    # Pseudobulk samples with too few cells
    low_cells <- which(agg$ncells < min_cells)
    if (verbose) message(.bluem("[DE] "), "Samples with insufficient cell numbers: ", length(low_cells))
    if (length(low_cells) > 0) agg <- agg[, -low_cells]

    # Labels with too few replicates in the conditions of the contrast,
    # and samples from other conditions in groups with too few replicates
    rep_tab <- table(agg[[labels]], agg[[condition]])
    low_rep_labels <- rownames(rep_tab)[apply(rep_tab[, conditions_used, drop = FALSE] < min_rep, 1, any)]
    group_size <- rep_tab[cbind(as.character(agg[[labels]]), as.character(agg[[condition]]))]
    drop_labels <- as.character(agg[[labels]]) %in% low_rep_labels
    drop_groups <- !agg[[condition]] %in% conditions_used & group_size < min_rep

    if (verbose && length(low_rep_labels) > 0) {
        message(.bluem("[DE] "),
                "The following labels were removed due to low replicates/low cells: \n",
                paste0("\t", paste(low_rep_labels, collapse = ", ")))
    }
    if (verbose && any(drop_groups & !drop_labels)) {
        message(.bluem("[DE] "),
                "Samples removed from conditions outside the contrast due to low replicates: ",
                sum(drop_groups & !drop_labels))
    }
    agg <- agg[, !(drop_labels | drop_groups)]

    missing_conditions <- setdiff(conditions_used, as.character(agg[[condition]]))
    if (length(missing_conditions) > 0) {
        stop(paste0(ep, "No pseudobulk samples are left for the condition(s): ",
                    paste(missing_conditions, collapse = ", "),
                    " after filtering. Consider lowering min_cells or min_rep."))
    }

    if (verbose) message(.bluem("[DE] "), "Performing differential expression analysis.")

    # Remove conditions that are missing in some labels and are not in the contrast
    if (!subset_conditions) {
        check_full_rank <- lapply(unique(as.character(agg[[labels]])), function(l) {
            curdata <- colData(agg[, agg[[labels]] == l])
            curdesign <- model.matrix(design, data = curdata)
            names(which(colSums(curdesign) == 0))
        })
        check_full_rank <- Reduce(union, check_full_rank)
        check_full_rank <- check_full_rank[startsWith(check_full_rank, condition)]
        remove_conditions <- setdiff(check_full_rank, paste0(condition, conditions_used))
        remove_conditions <- substring(remove_conditions, nchar(condition) + 1)
        agg <- agg[, !agg[[condition]] %in% remove_conditions]
    }
    colData(agg)[[condition]] <- droplevels(agg[[condition]])
    colData(agg)[[labels]] <- droplevels(agg[[labels]])

    dge <- pseudoBulkDGE(agg,
                         label = agg[[labels]],
                         condition = agg[[condition]],
                         col.data = colData(agg),
                         design = design,
                         contrast = contrast)

    # Report labels that could not be tested
    failed <- vapply(dge, function(x) isTRUE(metadata(x)$failed), logical(1))
    not_returned <- setdiff(levels(agg[[labels]]), names(dge))
    if (verbose && any(failed)) {
        message(.bluem("[DE] "),
                "The following labels could not be tested (results are NA): \n",
                paste0("\t", paste(names(dge)[failed], collapse = ", ")))
    }
    if (verbose && length(not_returned) > 0) {
        message(.bluem("[DE] "),
                "No results were returned for the following labels: \n",
                paste0("\t", paste(not_returned, collapse = ", ")))
    }

    for (i in seq_along(dge)) dge[[i]]$label <- names(dge)[i]

    dge <- lapply(dge, as.data.frame)

    attributes(dge)$conditions <- conditions_used
    attributes(dge)$inputs <- c(replicates, labels, condition)

    if (verbose) message(.bluem("[DE] "), "Done.")

    dge
}

#' Per-label GSEA
#'
#' Perform GSEA on a list of per-label DE results
#' 
#' @param dge a list of data frames with the differential expression results
#'    from the \code{doPBDGE} function.
#' @param pathways a named list of genesets/pathways to use in the GSEA
#' @param alpha numeric, the FDR threshold to use for filtering
#' 
#' @return a named list with GSEA results (\code{gsea}), and a matrix (\code{nesmat}) 
#' 	   where each column represents a label DE for which GSEA was performed, and each 
#'     row a geneset/pathway that is significant in at least one of the labels.
#' 
#' 
#' @export

doLabelGSEA <- function(dge, 
					    pathways, 
					    alpha = 0.05) {

	### Sanity checks
		# error prefix
		ep <- .redm("{cellula::GSEAMatrix()} - ")
  
	dependencies = data.frame("package" = c("fgsea"),
                              "repo" = c("BioC"))
    if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))

	if(is.null(names(pathways))) stop(paste0(ep, "pathways must be a named list"))
	if(!is(dge, "list")) stop(paste0(ep, "dge must be a list"))

	gsea_list = lapply(dge, function(x) {
		x = x[!is.na(x$logFC),]
		stats = x$logFC
		names(stats) = rownames(x)
		message("GSEA for ", unique(x$label), " using ", length(stats), " genes...\n")
		fgsea::fgsea(pathways = pathways, stats = stats)
	})
	
	names(gsea_list) = names(dge)
	nesmat_ok = makeNESmatrix(gsea_list, alpha = alpha)
	
	list("gsea" = gsea_list, 
		"nesmat" = nesmat_ok)
}

#' GSEA NES matrix
#'
#' Create a NES matrix from a list of GSEA results
#' 
#' @param gsea_list a list of GSEA results from \code{link{doLabelGSEA}}.
#' @param alpha numeric, the FDR threshold to use for filtering
#' 
#' @return a matrix where each column represents a label DE 
#'     for which GSEA was performed, and each row a geneset/pathway 
#'     that is significant in at least one of the labels.
#' 
#' @export

makeNESmatrix <- function(gsea_list, alpha = 0.05) {
  
  sig_path_union = Reduce(union, lapply(gsea_list, function(x) 
    x$pathway[which(x$padj < alpha)]))
  
  nesmat = matrix(nrow = length(sig_path_union), 
                  ncol = length(gsea_list), 
                  dimnames = list(sig_path_union, names(gsea_list)))
  
  for(i in seq_along(gsea_list)) {
    rownames(gsea_list[[i]]) = gsea_list[[i]]$pathway
    vec = rep(0, length(sig_path_union))
    names(vec) = sig_path_union
    vec[rownames(gsea_list[[i]])[which(gsea_list[[i]]$padj < alpha)]] = gsea_list[[i]]$NES[which(gsea_list[[i]]$padj < alpha)]
    nesmat[,names(gsea_list)[i]] = vec
  }
  
  nesmat[,apply(nesmat, 2, function(x) !all(is.na(x)))]
}

#' Per-label GSEA
#'
#' Perform GSEA on a list of per-label DE results
#' 
#' @param dge a list of data frames with the differential expression results
#'    from the \code{doPBDGE} function.
#' @param pathways a named list of genesets/pathways to use in the GSEA
#' @param alpha numeric, the FDR threshold to use for filtering
#' 
#' @return a named list with GSEA results (\code{gsea}), and a matrix (\code{nesmat}) 
#' 	   where each column represents a label DE for which GSEA was performed, and each 
#'     row a geneset/pathway that is significant in at least one of the labels.
#' 
#' 
#' @export

doLabelGSEA <- function(dge, 
					    pathways, 
					    alpha = 0.05) {

	### Sanity checks
		# error prefix
		ep <- .redm("{cellula::GSEAMatrix()} - ")
  
	dependencies = data.frame("package" = c("fgsea"),
                              "repo" = c("BioC"))
    if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))

	if(is.null(names(pathways))) stop(paste0(ep, "pathways must be a named list"))
	if(!is(dge, "list")) stop(paste0(ep, "dge must be a list"))

	gsea_list = lapply(dge, function(x) {
		x = x[!is.na(x$logFC),]
		stats = x$logFC
		names(stats) = rownames(x)
		message("GSEA for ", unique(x$label), " using ", length(stats), " genes...\n")
		fgsea::fgsea(pathways = pathways, stats = stats)
	})
	
	names(gsea_list) = names(dge)
	nesmat_ok = makeNESmatrix(gsea_list, alpha = alpha)
	
	list("gsea" = gsea_list, 
		"nesmat" = nesmat_ok)
}

#' GSEA NES matrix
#'
#' Create a NES matrix from a list of GSEA results
#' 
#' @param gsea_list a list of GSEA results from \code{link{doLabelGSEA}}.
#' @param alpha numeric, the FDR threshold to use for filtering
#' 
#' @return a matrix where each column represents a label DE 
#'     for which GSEA was performed, and each row a geneset/pathway 
#'     that is significant in at least one of the labels.
#' 
#' @export

makeNESmatrix <- function(gsea_list, alpha = 0.05) {
  
  sig_path_union = Reduce(union, lapply(gsea_list, function(x) 
    x$pathway[which(x$padj < alpha)]))
  
  nesmat = matrix(nrow = length(sig_path_union), 
                  ncol = length(gsea_list), 
                  dimnames = list(sig_path_union, names(gsea_list)))
  
  for(i in seq_along(gsea_list)) {
    rownames(gsea_list[[i]]) = gsea_list[[i]]$pathway
    vec = rep(0, length(sig_path_union))
    names(vec) = sig_path_union
    vec[rownames(gsea_list[[i]])[which(gsea_list[[i]]$padj < alpha)]] = gsea_list[[i]]$NES[which(gsea_list[[i]]$padj < alpha)]
    nesmat[,names(gsea_list)[i]] = vec
  }
  
  nesmat[,apply(nesmat, 2, function(x) !all(is.na(x)))]
}

