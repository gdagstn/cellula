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

			if(is.null(design)){
					design = formula(paste0("~0 + ", condition))
			}
			agg$label = colData(agg)[,labels]

			if(verbose) message(.bluem("[DE] "),"Performing differential expression analysis.")

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


#' Differential gene expression stripchart
#' 
#' Plots a stripchart of the differential gene expression results
#' 
#' @param dge a data frame with the differential expression results
#'    from the \code{doPBDGE} function
#' @param alpha numeric, the FDR threshold to use for filtering
#' @param lfc numeric, the log-fold change threshold to use for filtering
#' @param title character, the title of the plot
#' 
#' @return a stripchart with per-cluster DE results
#' 
#' @importFrom ggplot2 ggplot geom_vline geom_label ggtitle
#' @importFrom ggplot2 scale_color_gradient theme_bw .data aes
#' @importFrom ggbeeswarm geom_quasirandom
#' 
#' @export

plotDGEStripchart <- function(dge, 
							  alpha = 0.05, 
							  lfc = 1, 
							  title = NULL) {

		conds = attr(dge, "conditions")

		dge = do.call(rbind, dge)
		dge = dge[!is.na(dge$logFC),]

		ndegs = unlist(lapply(split(dge, dge$label), function(x) length(which(abs(x$logFC) > lfc & x$FDR < alpha))))
		ndegs = data.frame(label = names(ndegs), ndegs = ndegs)

		dge$label = factor(dge$label, levels = ndegs$label[order(ndegs$ndegs, decreasing = FALSE)])

		p = ggplot(dge[dge$FDR > alpha,], aes(x = .data[["logFC"]], y = .data[["label"]])) +
			geom_quasirandom(method = "tukey", 
							 color="gray", 
							 alpha = 0.5, 
							 orientation="y") +
			geom_quasirandom(data = dge[dge$FDR < alpha & abs(dge$logFC) > lfc,], 
							aes(x = .data[["logFC"]], color = -log10(.data[["FDR"]]), y = .data[["label"]]), 
							pch = 16, 
							inherit.aes = FALSE, 
							method = "tukey", 
							orientation = "y") +
			scale_color_gradient(low = "gray60", high = "orange") +
			geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.5) +
			geom_label(data = ndegs, mapping = aes(y = .data[["label"]], label = .data[["ndegs"]]), x = 0) + 
			theme_bw() 

		if(is.null(title)) {
			title = paste0(conds[1], " vs ", conds[2], " - log2(FC) distribution per label")
		}

	p + ggtitle(title)
}

#' Differential gene expression heatmap
#' 
#' Draws a heatmap of aggregated genes together with their DE results
#' 
#' @param sce a SingleCellExperiment object
#' @param genes character, the genes to include in the heatmap.
#'     These will be looked up in the rownames of the object, which are also the 
#' 	   rownames of the data frames in the DE results. At least 2 genes must be provided.
#' @param dge a data frame with the differential expression results from \code{\link{doPBDGE}}
#' @param exprs character, the expression values to use. Default is \code{"logcounts"}
#' @param alpha_include numeric, the FDR threshold to include a gene in the heatmap
#' 		Default is 1 (all genes included)
#' @param alpha_dot numeric, the FDR threshold to mark a gene as significant in
#' 		the heatmap using a dot. Default is 0.05
#' @param condition_pal a named vector with the colors for the conditions. Default is NULL
#' @param cluster_pal a named vector with the colors for the clusters. Default is NULL
#' @param title character, the title of the heatmap. Default is NA, no title will be displayed.
#' @param cluster_cols logical, should the columns be clustered? Default is \code{TRUE}
#' @param gaps logical, should gaps be included between the columns? Default is \code{FALSE}
#' @param scale logical, should the expression data be scaled by row? Default is \code{TRUE}
#' 
#' @return a composite heatmap of the pseudo-bulk aggregated genes with their DE results
#' 
#' @details this function generates a heatmap of genes using pseudo-bulk aggregated profiles
#'    taken by aggregating the \code{sce} object by the variables stored in the attributes
#'    of the \code{dge} object, and taking the assay specified by \code{exprs}. 
#'    Then, a second heatmap is drawn on the right hand side of the first heatmap, showing
#'    the DE results for the contrast stored in the attributes of the \code{dge} object.
#' 		
#'    The genes are filtered by the \code{alpha_include} threshold, and significant genes
#'    are labelled in the right hand side heatmap using a white dot if their FDR is below
#' 	  \code{alpha_dot}.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom scuttle aggregateAcrossCells
#' 
#' 
#' @export

plotDEHeatmap <- function(sce,
                          genes, 
                          dge, 
                          exprs = "logcounts", 
                          alpha_include = 1, 
                          alpha_dot = 0.05,
						  condition_pal = NULL,
						  cluster_pal = NULL,
                          title = NA,
                          cluster_cols = FALSE,
                          gaps = TRUE,
                          scale = TRUE) {

	### Sanity checks
		# error prefix
			ep = .redm("cellula::drawDEHeatmap() - ")
	if(!is(sce, "SingleCellExperiment")) 
		stop(paste0(ep, "Must provide a SingleCellExperiment object."))
	if(!is(genes, "character")) 
		stop(paste0(ep, "genes must be a character vector"))

  	dependencies = data.frame("package" = c("ComplexHeatmap", "circlize"),
                              	  "repo" = c("BioC", "CRAN"))
    if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))

  selected = intersect(genes, rownames(sce))
  
  	if(length(selected) <= 1) 
  		stop(paste0(ep,"Need at least 2 genes that are present in the object."))
  
  selected_ok = Reduce(union, lapply(dge, 
                                     function(x) {
                                       rownames(x)[rownames(x) %in% selected & x$FDR < alpha_include]
                                     })
  )
  
  selected_ok = selected_ok[!is.na(selected_ok)]    
  
  delist = lapply(names(dge), function(x) 
    data.frame(gene = selected_ok,                    
               lfc = dge[[x]][selected_ok, "logFC"],
               fdr = dge[[x]][selected_ok, "FDR"],
               cluster = x)
  )

  conditions_used = attr(dge, "conditions")

  conditions_vector = colData(sce)[,attr(dge, "inputs")[3]]
  replicates_vector = colData(sce)[,attr(dge, "inputs")[1]][which(conditions_vector %in% conditions_used)]
  labels_vector = colData(sce)[,attr(dge, "inputs")[2]][which(conditions_vector %in% conditions_used)]
  conditions_vector = conditions_vector[which(conditions_vector %in% conditions_used)]

  agg = aggregateAcrossCells(sce[,which(colData(sce)[,attr(dge, "inputs")[3]] %in% conditions_used)],
  							 ids = paste0(conditions_vector, "__",
										  replicates_vector, "__",
										  labels_vector),
									statistics = "sum", 
									use.assay.type = "counts")
  agg = computeLibraryFactors(agg)
  agg = logNormCounts(agg)

  agg$condition =  as.character(colData(agg)[,attr(dge, "inputs")[3]])
  agg$label = as.character(colData(agg)[,attr(dge, "inputs")[2]])
  agg$replicates =  as.character(colData(agg)[,attr(dge, "inputs")[1]])

  names(delist) = names(dge)
  dedf = do.call(rbind, delist)
  
  dat = assay(agg, exprs)[selected_ok,]
  dat = dat[,order(agg$label, agg$condition)]
  if(scale) dat = t(scale(t(dat)))
  
  res_lfc = lapply(delist, function(x) x$lfc)
  res_FDR =  lapply(delist, function(x) x$fdr)
  names(res_lfc) = names(res_FDR) = names(delist)
  
  pch_sigs = lapply(res_FDR, function(x) ifelse(x < alpha_dot, 21, NA))
  
  cr = max(abs(range(unlist(res_lfc), na.rm = TRUE)))
  color_range = c(-cr, 0, cr)
  
  pal_cond = condition_pal

  if(is.null(pal_cond)) {
	conds = attr(dge, "conditions")
	pal_cond = c("orange", "purple")
	names(pal_cond) = conds
  }

  pal_cluster = cluster_pal
	
	if(is.null(cluster_pal)) {
	pal_cluster = cellula:::.cpal_qual(length(dge))
	names(pal_cluster) = unique(names(dge))
	}

	cl_labs = agg$label[order(agg$label, agg$condition)]
	con_labs = agg$condition[order(agg$label, agg$condition)]

	if(gaps) colgaps = cl_labs else colgaps = NULL

	column_ha = ComplexHeatmap::HeatmapAnnotation(label = cl_labs,
												condition =  con_labs,
												col = list(label = pal_cluster,
															condition = pal_cond))
	cf = circlize::colorRamp2(color_range, 
							c(pal_cond[2], "gray", pal_cond[1]))
	
	ras = lapply(names(delist), function(x) {
		res_lfc[[x]][is.na(res_lfc[[x]])] = 0
		ComplexHeatmap::anno_simple(res_lfc[[x]], 
						which = "row",
						col = cf,
						pch = pch_sigs[[x]],
						pt_gp = grid::gpar(bg = "white"),
						pt_size = unit(1.5, "mm"))
	})
	names(ras) = names(delist)
	sig_anno = do.call(eval(parse(text="ComplexHeatmap::rowAnnotation")), ras)
	
	lgd = ComplexHeatmap::Legend(col_fun = cf, title = "log2(FC)")

	H = ComplexHeatmap::Heatmap(dat, 
								row_names_gp = grid::gpar(cex = 1), 
								column_names_gp = grid::gpar(cex = 1), 
								top_annotation = column_ha,
								cluster_columns = cluster_cols,
								right_annotation = sig_anno,
								column_split = colgaps,
								column_title = title,
								show_column_names = FALSE,
								name = "Z score"
					)
	
	ComplexHeatmap::draw(H, annotation_legend_list = list(lgd))
}

#' Per-label PCA plot
#' 
#' Plot a PCA plot of the aggregated data per label using
#' the first 2 components and top 1000 HVGs
#' 
#' @param sce a SingleCellExperiment object
#' @param replicates character, column in \code{colData(sce)}
#'     that contains the replicate information
#' @param labels character, column in \code{colData(sce)}
#'     that contains the label information
#' @param condition character, column in \code{colData(sce)}
#'     that contains the condition information
#' 
#' @return a PC1 vs PC2 faceted plot for pseudobulk profiles 
#'     where each panel corresponds to a label
#' 
#' @importFrom scuttle aggregateAcrossCells computeLibraryFactors logNormCounts 
#' @importFrom scran modelGeneVar getTopHVGs
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme_bw
#' @importFrom stats prcomp
#' 
#' @export 

plotLabelPCA <- function(sce, 
						replicates,
						labels,
						condition){

			conditions_vector = colData(sce)[,condition]
			replicates_vector = colData(sce)[,replicates]
			labels_vector = colData(sce)[,labels]

			agg = aggregateAcrossCells(sce,
									ids = paste0(conditions_vector, "__",
												 replicates_vector, "__",
												 labels_vector),
									statistics = "sum", 
									use.assay.type = "counts")

			agg = computeLibraryFactors(agg)
			agg = logNormCounts(agg)

			agg$label = colData(agg)[,labels]
			agg$condition = colData(agg)[,condition]
			agg$replicate = colData(agg)[,replicates]
			
			vargenes = modelGeneVar(agg, block = agg$label)

			hvglist = lapply(vargenes$per.block, getTopHVGs, n = 1000)
			names(hvglist) = names(vargenes$per.block)

			pcalist = lapply(names(hvglist), function(x) {
				prcomp(t(assay(agg, "logcounts")[hvglist[[x]], agg$label == x]), scale = TRUE)
			})

			pcalist = lapply(pcalist, function(x) as.data.frame(x$x[,1:2]))
			names(pcalist) = names(hvglist)
			for(i in seq_along(pcalist)) {
				pcalist[[i]]$label = names(pcalist)[i]
				pcalist[[i]]$condition = agg$condition[agg$label == names(pcalist)[i]]
				pcalist[[i]]$replicate = agg$replicate[agg$label == names(pcalist)[i]]
			}
			pcadf = do.call(rbind, pcalist)

	 ggplot(pcadf, aes(x = .data[["PC1"]], y = .data[["PC2"]])) +
			geom_point(aes(color = .data[["condition"]])) +
			facet_wrap("label")+
			theme_bw()
}

plotLabelMD <- function(dge, 
						ntop = 5, 
						alpha = 0.05, 
						lfc = 1) {

	for(i in seq_along(dge)) {
		dge[[i]] = dge[[i]][!is.na(dge[[i]]$logFC),]
		dge[[i]]$gene = rownames(dge[[i]])
	}

	top_genes = lapply(dge, function(x) {
		x = x[abs(x$logFC) > lfc & x$FDR < alpha,]
		x = x[order(abs(x$logFC), decreasing = TRUE),]
		x = x[seq_len(min(ntop, nrow(x))),]
		x
	})
	top_genes = do.call(rbind, top_genes)

	des = do.call(rbind, dge)
	des = des[order(des$FDR, decreasing = TRUE),]

	
	 ggplot(des, aes(x = .data[["logCPM"]], y = .data[["logFC"]])) +
			geom_point(aes(color = -log10(.data[["FDR"]]))) +
			geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2) +
			geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2) +
			geom_text_repel(data = top_genes, 
							aes(label = .data[["gene"]]), 
							box.padding = 0.5,
							size = 2.4,
							segment.size = 0.2) +
			scale_color_gradient(low = "gray60", high = "orange") +
			facet_wrap("label")+
			theme_bw()
}


# GSEA

GSEAmatrix <- function(dge, pathways, alpha = 0.05) {

  gsea_list = list()
  des = split(dge, dge$label)
  for(i in seq_along(des)) {
    des[[i]] = des[[i]][!is.na(des[[i]]$logFC),]
    stats = des[[i]]$logFC
    names(stats) = rownames(des[[i]])
    message("GSEA for ", names(des)[i], " using ", length(stats), " genes...\n")
    gsea_list[[i]] = fgsea::fgsea(pathways = pathways, stats = stats, )
  }
  names(gsea_list) = names(des)
  nesmat_ok = makeNESmatrix(gsea_list, alpha = alpha)
  list("gsea" = gsea_list, "nesmat" = nesmat_ok)
}

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
