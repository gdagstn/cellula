#' Pseudo-bulk differential expression analysis
#' 
#' @description This function performs differential expression
#'     analysis on a single-cell experiment with different replicates, 
#' 	   aggregated into pseudobulk on a per-cluster basis.
#' 
#' @param sce a SingleCellExperiment object
#' @param replicates character, column in \code{colData(sce)}
#'     that contains the replicate information
#' @param labels character, column in \code{colData(sce)}
#'     that contains the label information
#' @param condition character, column in \code{colData(sce)}
#'     that contains the condition information
#' @param contrast character, the contrast to be tested.
#'     Default is the second level of the condition minus the first level.
#'     See Details for how it is formatted.
#' @param design a design matrix for the differential expression analysis.
#' 	   Default is NULL, meaning the design matrix will be created automatically
#' 	   with an intercept term set to 0, e.g. \code{~0 + condition}.
#' @param verbose logical, display messages on progress? 
#'     Default is \code{TRUE}.
#' @param parallel_param a \code{BiocParallel} object specifying
#'     the parallelization backend
#' 
#' @return a data frame with the differential expression results
#'     through the \code{\link[scran]{pseudoBulkDGE}} function from \code{{scran}}.
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
#'    The contrast to be tested can be specified as well. If not, the default
#'    is the second level of the condition minus the first level. The format
#'    of the contrast is the one accepted by \code{\link[scran]{pseudoBulkDGE}}, i.e.
#'    \code{"conditionLevel1-conditionLevel2"} where condition is the name of
#'    the column containing the condition, and Level1 and Level2 are the two
#'    levels to compare.
#' 
#'    Default design is a formula with an intercept term set to 0, e.g.
#'    \code{~0 + condition} as it forces the user to explicitly declare in the
#'    contrast which is the reference level. However, if design is provided,
#'    the user can specify a design without intercept (e.g. \code{~condition})
#'    which means that the reference level for the factor will be taken as the
#'    reference level for the contrast too. Not knowing what the reference
#'    value for the factor is can cause the function to fail if the reference
#'    level is used in the \code{contrast} argument.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom scran pseudoBulkDGE
#' @importFrom scuttle computeLibraryFactors logNormCounts aggregateAcrossCells
#' @importFrom BiocParallel SerialParam
#' @importFrom stats formula
#' 
#' @export

doPBDGE <- function(sce, 
					replicates,
					labels,
					condition,
					contrast = NULL,
					design = NULL,
					verbose = TRUE,
					parallel_param = SerialParam()) {

### Sanity checks

			ep = .redm("cellula::doPBDGE() - ")

			if(!is(sce, "SingleCellExperiment")) 
			stop(paste0(ep, "Must provide a SingleCellExperiment object."))

			if(!is.character(replicates)) 
			stop(paste0(ep, "replicates must be a character."))
			if(!is.character(labels)) 
			stop(paste0(ep, "labels must be a character."))
			if(!is.character(condition)) 
			stop(paste0(ep, "condition must be a character."))

			if(!condition %in% colnames(colData(sce))) 
			stop(paste0(ep, "The condition column was not found in the colData."))
			if(!replicates %in% colnames(colData(sce))) 
			stop(paste0(ep, "The replicates column was not found in the colData."))
			if(!labels %in% colnames(colData(sce))) 
			stop(paste0(ep, "The labels column was not found in the colData."))

			if(!is.null(contrast) && !is.character(contrast)) 
			stop(paste0(ep, "The contrast must be a character."))
			if(!is.null(design) && !inherits(design, "formula")) 
			stop(paste0(ep, "The design must be a formula."))
	
			if(verbose) message(.bluem("[DE] "),"Aggregating into pseudobulk.")

			conditions_vector = colData(sce)[,condition]
			replicates_vector = colData(sce)[,replicates]
			labels_vector = colData(sce)[,labels]

			agg = aggregateAcrossCells(sce,
									ids = paste0(conditions_vector, "__",
												 replicates_vector, "__",
												 labels_vector),
									statistics = "sum", 
									use.assay.type = "counts",
									BPPARAM = parallel_param)

			agg$condition = factor(colData(agg)[,condition])

			if(is.null(contrast)) {
				message(.bluem("[DE] "),"Contrast was not specified. Defaulting to: ", 
				condition, levels(agg$condition)[2], "-", condition, levels(agg$condition)[1])

				contrast = paste0(condition, levels(agg$condition)[2], "-", condition, levels(agg$condition)[1])
			}
			agg$label = colData(agg)[,labels]

			if(verbose) message(.bluem("[DE] "),"Performing differential expression analysis.")

			if(is.null(design)) design = formula(paste0("~0 + ", condition))

			dge = pseudoBulkDGE(agg, 
								label = agg$label, 
								condition = agg$condition,
								col.data = colData(agg), 
								design = design,
								contrast = contrast)

			if(verbose) message(.bluem("[DE] "),"Adding mean log-normalized counts per gene.")

			agg_mean = 	aggregateAcrossCells(sce,
									ids = paste0(labels_vector),
									statistics = "mean", 
									use.assay.type = "counts",
									BPPARAM = parallel_param)

			agg_mean = computeLibraryFactors(agg_mean)
			agg_mean = logNormCounts(agg_mean)

			for(i in names(dge)) {
				dge[[i]]$mean = assay(agg_mean, "logcounts")[,i]
			}

			if(verbose) message(.bluem("[DE] "),"Done.")

			for(i in seq_along(dge)) dge[[i]]$label = names(dge)[i]

			#des = do.call(rbind, dge)

			des
}
