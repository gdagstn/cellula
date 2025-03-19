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
#' @param min_cells numeric, minimum number of cells per condition
#'     to keep the pseudobulk sample. Default is 10.
#' @param subset_conditions logical, should the object be subset to contain 
#'	   only the conditions specified in the contrast? 
#'     Default is \code{FALSE}.
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
#' 	  An option is given to subset the object to only contain the two conditions
#'    included in the contrast. The difference between subsetting and not 
#' 	  lies in the fact that dispersion estimates will only be calculated for
#'    the samples included in the subset object, and not in the whole object.
#'    This can result in different results (e.g. fewer or more DE genes). 
#'    According to current best practices it is NOT recommended to subset
#'    the object, but the option is given.
#' 
#'    The contrast to be tested can be specified as well. If not, the default
#'    is the second level of the condition minus the first level. The format
#'    of the contrast is the one accepted by \code{\link[scran]{pseudoBulkDGE}}, i.e.
#'    \code{"conditionLevel1-conditionLevel2"} where condition is the name of
#'    the column containing the condition, and Level1 and Level2 are the two
#'    levels to compare.
#' 
#'    The default design is a formula with an intercept term set to 0, e.g.
#'    \code{~0 + condition} as it forces the user to explicitly declare in the
#'    contrast which is the reference level. However, if design is provided,
#'    the user can specify a design without intercept (e.g. \code{~condition})
#'    which means that the reference level for the factor will be taken as the
#'    reference level for the contrast too. Not knowing what the reference
#'    value for the factor is can cause the function to fail if the reference
#'    level is used in the \code{contrast} argument.
#' 
#' @importFrom SummarizedExperiment colData 
#' @importFrom SingleCellExperiment sizeFactors
#' @importFrom scran pseudoBulkDGE 
#' @importFrom scuttle computeLibraryFactors logNormCounts aggregateAcrossCells
#' @importFrom scuttle librarySizeFactors
#' @importFrom BiocParallel SerialParam
#' @importFrom stats formula
#' 
#' @export

doPBDGE <- function(sce, 
					replicates,
					labels,
					condition,
					min_cells = 10,
					subset_conditions = FALSE,
					contrast = NULL,
					design = NULL,
					verbose = TRUE,
					parallel_param = SerialParam()) {

	### Sanity checks
		# error prefix
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

			conditions_used = unlist(strsplit(gsub(condition, "", contrast), "-"))
			if(length(conditions_used) == 1) conditions_used = c(conditions_used, levels(sce$condition)[1])
			
			if(subset_conditions) {
				if(verbose) message(.bluem("[DE] "),
							"The object will be subset to contain the following conditions only: \n", 
							paste0("\t", paste(conditions_used, collapse = ", ")))

				sce = sce[,sce$condition %in% conditions_used]
			}

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
			agg$label = factor(colData(agg)[,labels])

			# Filtering
			remove_ids = agg$ids[which(agg$ncells < min_cells & agg$condition %in% conditions_used)]

			if(verbose) message(.bluem("[DE] "), "Samples with insufficient cell numbers: ", length(remove_ids))

			all_tab = as.data.frame(table(as.data.frame(colData(agg))[,c("condition", "label")]))
			all_tab_current = all_tab[all_tab$condition %in% conditions_used,]
			remove_tab = as.data.frame(table(as.data.frame(colData(agg))[remove_ids, c("condition", "label")]))

			remove_lowrep = Reduce(union, lapply(split(all_tab_current, all_tab_current$condition), 
												function(x) as.character(x$label[which(x$Freq < 3)])))

			remove_tab$Freq[remove_tab$label %in% remove_lowrep] = Inf

			remove_tab_ok = all_tab[which((all_tab$Freq - remove_tab$Freq) < 3),]

			remove_samples = unlist(apply(remove_tab_ok, 1, function(x) 
				which(agg$condition == x[1] & agg$label == x[2])))

			before = unique(agg$label)
			if(length(remove_samples) > 0) {
				agg = agg[, -1 * remove_samples]
				difference = setdiff(before, unique(agg$label))
				if(verbose) message(.bluem("[DE] "),
							"The following labels were removed due to low replicates/low cells: \n", 
							paste0("\t", paste(difference, collapse = ", ")))
			}

			if(verbose) message(.bluem("[DE] "),"Performing differential expression analysis.")
			
			if(is.null(design)){
					design = formula(paste0("~0 + ", condition))
			}

			dge = pseudoBulkDGE(agg, 
								label = agg$label, 
								condition = agg$condition,
								col.data = colData(agg), 
								design = design,
								contrast = contrast)

			for(i in seq_along(dge)) dge[[i]]$label = names(dge)[i]
			
			dge = lapply(dge, as.data.frame)
			
			attributes(dge)$conditions = conditions_used
			attributes(dge)$inputs = c(replicates, labels, condition)

			if(verbose) message(.bluem("[DE] "),"Done.")

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

