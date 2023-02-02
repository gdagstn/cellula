#' Plot cluster modularity
#'
#' Plots a heatmap showing pairwise cluster modularity from a SCE object
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "modularity_SNN_100"
#' @param type character, one of "heatmap" or "graph" for different plot types
#'
#' @return a heatmap of pairwise modularity
#'
#' @importFrom pheatmap pheatmap
#' @importFrom colorspace sequential_hcl
#' @importFrom igraph graph_from_adjacency_matrix E layout_with_lgl
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export


plotModularity <- function(sce, name, type = "heatmap") {
  name = paste0("modularity_", name)
  modmat = log2(metadata(sce)[[name]] + 1)

  if(type == "heatmap"){
    pheatmap(mat = modmat,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = sequential_hcl(palette = "Sunset", n = 25),
             main = "log2(modularity ratio + 1)",
             border_color = NA
    )
  } else if(type == "graph") {

    cgr <- graph_from_adjacency_matrix(modmat,
                                       mode="upper",
                                       weighted=TRUE,
                                       diag=FALSE)

    set.seed(420)
    plot(cgr,
         edge.width=E(cgr)$weight*5,
         layout=igraph::layout_with_lgl,
         main = name)
  }
}


#' Plot approximate silhouette widths
#'
#' Plots a beesswarm plot showing the approximate width of each cell within
#' each cluster
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "silhouette_SNN_100"
#'
#' @return a beeswarm plot of silhouette widths
#'
#' @importFrom ggplot2 ggplot aes_string theme_bw
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details This plotting function has been directly taken from the OSCA book
#'     ("Orchestrating Single Cell Analysis with Bioconductor").
#'
#' @export

plotSilhouette <- function(sce, name) {
  name = paste0("silhouette_", name)
  ggplot(metadata(sce)[[name]], aes_string(x="cluster", y="width", colour="closest")) +
    geom_quasirandom(method="smiley") +
    theme_bw()
}


#' Plot UMAP
#'
#' Plots the UMAP from a SingleCellExperiment object
#'
#' @param sce a SingleCellExperiment object
#' @param umap_slot character, the name in the reducedDim slot of the SCE object, e.g.
#'     "UMAP_Harmony"
#' @param color_by character, column name in the colData slot of the SCE object, 
#'     e.g. "cluster". Will be used to assign colors. Automatically detects 
#'     whether the variable is categorical or continuous. 
#' @param shape_by character, column name in the colData slot of the SCE object,
#'     e.g. "cluster". Will be used to assign shapes. Can only be used for 
#'     categorical variables.
#' @param group_by character, column name in the colData slot of the SCE object,
#'     e.g. "individual". Will be used to facet the plot. Can only be used for
#'     categorical variables.
#' @param label_by character, column name in the colData slot of the SCE object,
#'     e.g. "cluster". Will be used to add labels to the plot. Can only be used
#'     for categorical variables.
#' @param color_palette a character string containing colors to be used. Default
#'     is NULL, meaning an automatic palette will be generated.
#' 
#' @return a ggplot object showing the UMAP coordinates colored, shaped and 
#'     faceted as set by the arguments. 
#'
#' @importFrom SingleCellExperiment reducedDimNames reducedDim 
#' @importFrom SummarizedExperiment colData
#' @importFrom colorspace sequential_hcl
#' @importFrom qualpalr qualpal
#' @importFrom ggplot2 ggplot aes scale_colour_gradientn guides geom_point
#' @importFrom ggplot2 guide_legend .data theme_void geom_segment geom_text coord_fixed 
#' @importFrom ggplot2 facet_wrap vars arrow unit geom_label scale_colour_manual theme
#' @importFrom rlang sym 
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export

plot_UMAP <- function(sce,  
                      umap_slot = "UMAP", 
                      color_by = NULL, 
                      shape_by = NULL, 
                      group_by = NULL, 
                      label_by = NULL,
                      color_palette = NULL) {
 
  ## Sanity checks
  # Error prefix
  ep = "{papplain::plot_UMAP()} - "
  
  # Checks
  if(!is(sce, "SingleCellExperiment")) 
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(!umap_slot %in% reducedDimNames(sce)) 
    stop(paste0(ep, "Could not find a slot named \"", umap_slot, "\" in the SingleCellExperiment object"))
  if(!is.null(color_by)){
    if(!color_by %in% colnames(colData(sce))) 
      stop(paste0(ep, color_by, " column not found in the colData of the SingleCellExperiment object"))
  }
  if(!is.null(shape_by)){
    if(!shape_by %in% colnames(colData(sce))) 
      stop(paste0(ep,shape_by, " column not found in the colData of the SingleCellExperiment object"))
  }
  if(!is.null(label_by)){
    if(!label_by %in% colnames(colData(sce))) 
      stop(paste0(label_by, " column not found in the colData of the SingleCellExperiment object"))
    if(!is(colData(sce)[,label_by], "character") & !is(colData(sce)[,label_by], "factor"))
      stop(paste0(ep, "Cannot label by a numeric value, convert it to factor first."))
  }
  if(!is.null(group_by)){
    if(!group_by %in% colnames(colData(sce))) 
      stop(paste0(ep, group_by, " column not found in the colData of the SingleCellExperiment object"))
    if(!is(colData(sce)[,group_by], "character") & !is(colData(sce)[,group_by], "factor"))
      stop(paste0(ep, "Cannot group by a numeric value, convert it to factor first."))
  }
  
  # Rescaling
  udf = as.data.frame(reducedDim(sce, umap_slot)[,1:2])
  udf[,1:2] = apply(udf[,1:2], 2, function(x) (x-min(x))/diff(range(x)))
  colnames(udf) = c("x", "y")

  # Check which aesthetics to include
  include = list(color_by, shape_by, group_by, label_by)
  include = unique(unlist(include[!is.null(include)]))
  
  for(i in include) udf[[i]] = colData(sce)[,i]
  
  colnames(udf)[seq_len(length(include)) + 2] = include
  
  classes = unlist(lapply(udf[,include], class))
  
  names(classes) = include
  aes_umap = aes(x = .data[["x"]], y = .data[["y"]])
  
  # Add color aesthetic if needed and decide what type of scale
  if(!is.null(color_by)) {
    aes_umap$colour = aes(colour = .data[[color_by]])$colour
    if(classes[color_by] == "numeric") {
      udf = udf[order(udf[,color_by]),]
      if(is.null(color_palette)) {
        pal = sequential_hcl(palette = "Sunset", n = 25)
      } else {
        pal = color_palette
      }
      cscale = scale_colour_gradientn(colours = pal) 
      cguides = NULL
    } else if(classes[color_by] == "character") {
      udf[,color_by] = factor(udf[,color_by])
      if(is.null(color_palette)) {
        pal = qualpal(n = length(unique(udf[,color_by])))$hex
      } else {
        pal = color_palette
      }
      cscale = scale_colour_manual(values = pal) 
      cguides = guides(color = guide_legend(override.aes = list(size=2)))
    } else if(classes[color_by] == "factor") {
      if(is.null(color_palette)) {
        pal = qualpal(n = length(unique(udf[,color_by])))$hex
      } else {
        pal = color_palette
      }
      cscale = scale_colour_manual(values = pal) 
      cguides = guides(color = guide_legend(override.aes = list(size=2)))
    }
  }
  
  # Add shape aesthetic if needed 
  if(!is.null(shape_by)) {
    if(classes[shape_by] == "numeric") {
      stop("Cannot map `shape` aesthetic to a numeric value")
    } else if(classes[shape_by] == "character"){
      udf[,shape_by] = factor(udf[,shape_by])
    }
    aes_umap$shape = aes(shape = .data[[shape_by]])$shape
  }
  
  if(!is.null(label_by)) {
    
    labels = unique(udf[,label_by])
    medoids = lapply(labels, function(x) {
      t(apply(udf[udf[,label_by] == x, c("x", "y")], 2, median))
    })
    label_df = as.data.frame(do.call(rbind, medoids))
    label_df[,label_by] = labels
    colnames(label_df) = c("x", "y", "label")
    
  }
  
  if(!is.null(group_by)) {
    if(is(colData(sce)[,group_by], "character"))
      colData(sce)[,group_by] = factor(colData(sce)[,group_by])
  }
  
  # Plot construction
  p = ggplot(udf, mapping = aes_umap) + 
    geom_point(size = 0.7, shape = 16) + 
    theme_void() + 
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) + 
    geom_segment(aes(x = 0, xend = 0.15, y = 0, yend = 0), 
                 color = "black",
                 linewidth = 0.3,
                 arrow = arrow(length=unit(0.1,"cm"), 
                               type = "closed")) + 
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0.15), 
                 color = "black",
                 linewidth = 0.3,
                 arrow = arrow(length=unit(0.1,"cm"), 
                               type = "closed")) + 
    geom_text(aes(x = -0.03, y = 0.02), 
              label = "UMAP 2", 
              angle = 90, 
              size = 2,
              hjust = "left",
              vjust = "top",
              color = "black") +
    geom_text(aes(x = 0.02, y = -0.03), 
              label = "UMAP 1", 
              size = 2,
              hjust = "left",
              vjust = "bottom",
              color = "black") +
    coord_fixed()
  
  if(!is.null(color_by)) {
    p = p + cscale + cguides
  }
  
  if(!is.null(label_by)) {
    p = p + geom_label(data = label_df, 
                       mapping = aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]]), 
                       inherit.aes = FALSE, size = 2)
  }
  
  if(!is.null(group_by)) {
    p = p + facet_wrap(vars(!!sym(group_by)))
  }
  
  return(p)
}

#' Plot a dot plot
#'
#' Plots a dot plot with gene expression from a SingleCellExperiment object 
#'
#' @param sce a SingleCellExperiment object
#' @param genes character string, the genes to be plotted (matched with rownames)
#' @param group_by character, column name in the colData slot of the SCE object, 
#'     e.g. "cluster". Will be used to assign calculate proportions. 
#'     Must be categorical (factor or coercible character). 
#' @param cluster_genes logical, should genes be clustered? Default is TRUE
#' @param cluster_groups logical, should groups be clustered? Default is TRUE
#' @param expres_use character, name of the `assay` slot in the SingleCellExperiment
#' @param color_palette a character string containing colors to be used. Default
#'     is NULL, meaning an automatic palette will be generated.
#' @param format character, one of "wide" (genes are columns, groups are rows) 
#'     or "tall" (genes are rows, grups are columns)
#' 
#' @return a ggplot object showing average gene expression and proportion of 
#'    expressing cells grouped by a grouping variable such as cluster.
#'
#' @importFrom SummarizedExperiment assay colData assayNames
#' @importFrom colorspace sequential_hcl
#' @importFrom stats dist hclust
#' @importFrom ggthemes theme_few
#' @importFrom ggplot2 ggplot aes scale_colour_gradientn geom_point 
#' @importFrom ggplot2 theme element_text labs element_blank element_line
#' @importFrom methods is
#'
#' @export

plot_dots <- function(sce, 
                     genes, 
                     group_by, 
                     cluster_genes = TRUE, 
                     cluster_groups = TRUE, 
                     exprs_use = "logcounts",
                     threshold = 0,
                     color_palette = NULL,
                     format = "wide") {
  
  
  ## Sanity checks
  # Error prefix
  ep = "{papplain::plot_Dots()} - "
  
  # Checks
  if(!is(sce, "SingleCellExperiment")) 
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(is.null(genes)) stop(paste0(ep, "Must provide genes"))
  if(!exprs_use %in% assayNames(sce)) 
    stop(paste0(ep, "Could not find an assay named \"", exprs_use, "\" in the SingleCellExperiment object"))
  if(is.null(group_by))
    stop(paste0(ep, "You must define a variable in group_by"))
  if(!group_by %in% colnames(colData(sce))) 
      stop(paste0(ep, group_by, " column not found in the colData of the SingleCellExperiment object"))
  if(!is(threshold, "numeric"))
      stop(paste0(ep, "Threshold must be numeric"))
  if(!is(colData(sce)[,group_by], "character") & !is(colData(sce)[,group_by], "factor"))
      stop(paste0(ep, "Cannot group by a numeric value, convert it to factor first."))
  if(any(!genes %in% rownames(sce))){
    found = intersect(genes, rownames(sce))
    if(length(found) == 0) {
      stop(paste0(ep, "None of the genes were found in the SingleCellExperiment object"))
    } else {
      message(paste0(ep, "Warning: some genes (", (length(genes) -length(found)), ") were not found in the dataset."))
      genes = found
    }
  }
  
  sce_byclust = lapply(unique(colData(sce)[,group_by]), function(x) {
    cur = assay(sce, exprs_use)[genes, which(colData(sce)[,group_by] == x)]
    props = apply(cur, 1, function(x) sum(x > threshold))/ncol(cur)
    aves = apply(cur, 1, function(x) mean(x))
    final_df = data.frame("proportion" = props, 
                          "mean_expression" = aves, 
                          "cluster" = x,
                          "gene" = genes)
    return(final_df)
  })
  
  if(cluster_genes) {
    scdf_props_genes = do.call(cbind, lapply(sce_byclust, function(x) x$proportion))
    scdf_exp_genes = do.call(cbind, lapply(sce_byclust, function(x) x$mean_expression))
    hc_props_genes = hclust(dist(scdf_props_genes))$order
    hc_exp_genes = hclust(dist(scdf_exp_genes))$order
  }
  
  if(cluster_groups) {
    scdf_props_clusters = do.call(rbind, lapply(sce_byclust, function(x) x$proportion))
    scdf_exp_clusters = do.call(rbind, lapply(sce_byclust, function(x) x$mean_expression))
    hc_props_clusters = hclust(dist(scdf_props_clusters))$order
    hc_exp_clusters = hclust(dist(scdf_exp_clusters))$order
  }
  
  scdf = do.call(rbind, sce_byclust)
  
  if(cluster_genes) scdf$gene = factor(scdf$gene, levels = genes[hc_exp_genes])
  if(cluster_groups) scdf$cluster = factor(scdf$cluster, levels = unique(colData(sce)[,group_by])[hc_exp_clusters])
  
  scdf$mean_expression[scdf$mean_expression == 0] <- NA
  scdf$proportion[scdf$proportion == 0] <- NA
  
  if(is.null(color_palette)) {
    pal = rev(sequential_hcl(n = 40, palette = "YlGnBu"))[10:40]
  } else {
    pal = color_palette
  }
  cscale = scale_colour_gradientn(colours = pal) 
  
  if(format == "wide") {
    p = ggplot(scdf, aes(y = cluster, x = gene, size = proportion, colour = mean_expression)) +
      geom_point(shape = 16, na.rm = TRUE) + 
      cscale + 
      theme_few() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            panel.border = element_blank(), 
            axis.line = element_line(colour = "black")) + 
      labs(colour = "Mean expression", size = "Proportion")
  } else if(format == "tall") {
    p = ggplot(scdf, aes(y = gene, x = cluster, size = proportion, colour = mean_expression)) +
      geom_point(shape = 16, na.rm = TRUE) + 
      cscale + 
      theme_few() + 
      theme(axis.text.y = element_text(face = "italic"),
            panel.border = element_blank(), 
            axis.line = element_line(colour = "black")) + 
      labs(colour = "Mean expression", size = "Proportion")
  }
  
  return(p)
}

