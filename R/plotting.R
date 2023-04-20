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
#' @importFrom ggplot2 ggplot aes_string theme_classic
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
    theme_classic()
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
#' @param point_size numeric, the size of the points in the plot. Default is 0.7.
#' @param label_size numeric, the size of the font for the labels. Default is 2.
#' @param color_palette a character string containing colors to be used. Default
#'     is NULL, meaning an automatic palette will be generated.
#' @param rescale logical, should coordinates be rescaled between 0 and 1? Default is TRUE.
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
#' @importFrom ggplot2 margin theme_minimal
#' @importFrom ggrepel geom_label_repel
#' @importFrom rlang sym
#' @importFrom stats median complete.cases
#' @importFrom methods is
#'
#' @export

plot_UMAP <- function(sce,
                      umap_slot = "UMAP",
                      color_by = NULL,
                      shape_by = NULL,
                      group_by = NULL,
                      label_by = NULL,
                      point_size = 0.7,
                      label_size = 2,
                      color_palette = NULL,
                      rescale = TRUE) {

  ## Sanity checks
  # Error prefix
  ep = "{cellula::plot_UMAP()} - "

  # Checks
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(!umap_slot %in% reducedDimNames(sce))
    stop(paste0(ep, "Could not find a slot named \"", umap_slot, "\" in the SingleCellExperiment object"))
  if(!is.null(color_by)){
    if(!color_by %in% colnames(colData(sce)) & !color_by %in% rowData(sce)$Symbol & !color_by %in% rowData(sce)$ID)
      stop(paste0(ep, color_by, " is not a column or a rowData element of the SingleCellExperiment object"))
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


  udf = as.data.frame(reducedDim(sce, umap_slot)[,1:2])
  udf = udf[complete.cases(udf),]

  # Rescaling
  if(rescale) {
    orig_range = lapply(udf, range)
    names(orig_range) = c("x", "y")
    udf[,1:2] = apply(udf[,1:2], 2, function(x) (x-min(x))/diff(range(x)))
  }

  colnames(udf) = c("x", "y")

  if(!is.null(color_by) & (color_by %in% rowData(sce)$Symbol | color_by %in% rowData(sce)$ID)) {
    feature = which(rowData(sce)$Symbol == color_by | rowData(sce)$ID == color_by)
    colData(sce)[,color_by] = as.numeric(assay(sce[feature,], "logcounts"))
  }

  # Check which mappings to include
  include = list(color_by, shape_by, group_by, label_by)
  include = unique(unlist(include[!is.null(include)]))

  for(i in include) udf[[i]] = colData(sce)[rownames(udf),i]

  colnames(udf)[seq_len(length(include)) + 2] = include

  classes = unlist(lapply(udf[,include,drop=FALSE], class))

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
    } else if(classes[color_by] %in% c("character", "logical")) {
      udf[,color_by] = factor(udf[,color_by])

      if(is.null(color_palette)) {
        if(length(levels(udf[,color_by])) == 1) {
          pal = "gray"
        } else {
          pal = qualpal(n = length(unique(udf[,color_by])))$hex
        }
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
    if(is(udf[,group_by], "character"))
      udf[,group_by] = factor(udf[,group_by])
  }

  # Plot construction
  arrow_1 = data.frame(x = 0, xend = 0.15, y = 0, yend = 0)
  arrow_2 = data.frame(x = 0, xend = 0, y = 0, yend = 0.15)
  text_1 = data.frame(x = -0.03, y = 0.02)
  text_2 = data.frame(x = 0.02, y = -0.03)

  p = ggplot(udf, mapping = aes_umap) +
    geom_point(size = point_size, shape = 16) +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    geom_segment(data = arrow_1,
                 mapping = aes(x = .data[["x"]],
                               y = .data[["y"]],
                               xend = .data[["xend"]],
                               yend = .data[["yend"]]),
                 inherit.aes = FALSE,
                 color = "black",
                 linewidth = 0.3,
                 arrow = arrow(length=unit(0.1,"cm"),
                               type = "closed")) +
    geom_segment(data = arrow_2,
                 mapping = aes(x = .data[["x"]],
                               y = .data[["y"]],
                               xend = .data[["xend"]],
                               yend = .data[["yend"]]),
                 color = "black",
                 inherit.aes = FALSE,
                 linewidth = 0.3,
                 arrow = arrow(length=unit(0.1,"cm"),
                               type = "closed")) +
    geom_text(data = text_1,
              aes(x = .data[["x"]],
                  y = .data[["y"]]),
              inherit.aes = FALSE,
              label = paste0(umap_slot, " 2"),
              angle = 90,
              size = 2,
              hjust = "left",
              vjust = "top",
              color = "black") +
    geom_text(data = text_2,
              aes(x = .data[["x"]],
                  y = .data[["y"]]),
              inherit.aes = FALSE,
              label = paste0(umap_slot, " 1"),
              size = 2,
              hjust = "left",
              vjust = "bottom",
              color = "black") +
    coord_fixed()

  if(!is.null(color_by)) {
    p = p + cscale + cguides
  }

  if(!is.null(label_by)) {
    p = p + geom_label_repel(data = label_df,
                             mapping = aes(x = .data[["x"]],
                                           y = .data[["y"]],
                                           label = .data[["label"]]),
                             inherit.aes = FALSE, size = label_size)
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
#' @param exprs_use character, name of the `assay` slot in the SingleCellExperiment
#' @param threshold numeric, expression threshold below which a gene is not expressed
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
  ep = "{cellula::plot_Dots()} - "

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

  if(cluster_genes) {
    scdf$gene = factor(scdf$gene, levels = genes[hc_exp_genes])
  } else scdf$gene = factor(scdf$gene, levels = rev(genes))
  if(cluster_groups) {
    scdf$cluster = factor(scdf$cluster, levels = unique(colData(sce)[,group_by])[hc_exp_clusters])
  }

  scdf$mean_expression[scdf$mean_expression == 0] <- NA
  scdf$proportion[scdf$proportion == 0] <- NA

  if(is.null(color_palette)) {
    pal = rev(sequential_hcl(n = 40, palette = "YlGnBu"))[10:40]
  } else {
    pal = color_palette
  }
  cscale = scale_colour_gradientn(colours = pal)

  if(format == "wide") {
    p = ggplot(scdf, aes(y = .data[["cluster"]],
                         x = .data[["gene"]],
                         size = .data[["proportion"]],
                         colour = .data[["mean_expression"]])
               ) +
      geom_point(shape = 16, na.rm = TRUE) +
      cscale +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
            ) +
      labs(colour = "Mean expression", size = "Proportion")
  } else if(format == "tall") {
    p = ggplot(scdf, aes(y = .data[["gene"]],
                         x = .data[["cluster"]],
                         size = .data[["proportion"]],
                         colour = .data[["mean_expression"]])
               ) +
      geom_point(shape = 16, na.rm = TRUE) +
      cscale +
      theme_minimal() +
      theme(axis.text.y = element_text(face = "italic"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
            ) +
      labs(colour = "Mean expression", size = "Proportion")
  }

  return(p)
}

#' @importFrom scater retrieveCellInfo

.makeColdataDF <- function(sce, columns) {
  do.call(cbind, lapply(columns, function(x) {
    ret = retrieveCellInfo(sce, by = x)
    col = data.frame(ret$value, row.names = colnames(sce))
    colnames(col) = ret$name
    return(col)
  })
  )
}

#' Plot cell metadata
#'
#' Plot data from the colData of the SCE object
#'
#' @param sce a SingleCellExperiment object
#' @param y character, the column from colData(sce) whose values need to be plotted
#' @param x character, the column from colData(sce) to be used as a grouping variable
#'     default is NULL, which means all points belong to one violin plot only
#' @param color_by character, the column name from colData(sce) to be used as a
#'     colouring variable. If NULL (default), violins (and points) will not be colored.
#' @param group_by character, the column name from colData(sce) to be used as a
#'     facetting variable. Default is NULL.
#' @param color_palette character string of colors for categorical plots or heatmap.
#'     Default is NULL.
#' @param contour logical, should contours be plotted on the scatterplot?
#'     default is TRUE.
#' @param clustered logical, should rows and columns be clustered in the confusion
#'     matrix heatmap? Default is TRUE.
#' @param plot_cells logical, should cells be plotted as well? Default is FALSE
#'     for speed.
#'
#' @returns a ggplot with different types of visualization depending on the
#'     classes of the columns chosen. If y is a numeric and x is NULL, or a
#'     categorical (character/facrtor), the function plots a violin + boxplot divided by
#'     x and/or by color_by. If y and x are both numeric, the function plots a
#'     scatterplot with optional 2D kernel contouring. If y and x are both
#'     categorical, the function plots a heatmap of the confusion matrix,
#'     showing pairwise Jaccard indices. All plots can be facetted.
#'
#' @importFrom ggplot2 facet_wrap
#' @importFrom rlang sym
#' @importFrom methods is
#'
#' @export

plot_Coldata <- function(sce,
                         y,
                         x = NULL,
                         color_by = NULL,
                         group_by = NULL,
                         color_palette = NULL,
                         contour = TRUE,
                         clustered = TRUE,
                         plot_cells = FALSE) {
  ## Sanity checks
  # Error prefix
  ep = "{cellula::plot_Coldata()} - "

  # Checks
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))

  # Prepare data
  include = c(x, y, color_by, group_by)
  include = unique(include[!is.null(include)])

  if(any(!include %in% colnames(colData(sce))))
    stop(paste0(ep, "Some columns were not found in the colData"))

  df <- .makeColdataDF(sce, include)

  if(!is.null(x) & is(df[,x], "numeric") & is(df[,y], "numeric")){
    p <- .makeScatter(df, x, y, color_by, color_palette, contour)
  } else if((!is(df[,x], "numeric") & !is.null(x)
             & is(df[,y], "numeric"))
            | (is.null(x) & is(df[,y], "numeric"))) {
    p <- .makeViolin(df, x, y, plot_cells, color_by, color_palette)
  } else if(!is.null(x)
            & (is(df[,x], "character") | is(df[,x], "factor"))
            & (is(df[,y], "character") | is(df[,y], "factor"))) {
    df[,x] = factor(df[,x])
    df[,y] = factor(df[,y])
    p <- .makeHeatmap(df, x, y, clustered, color_palette)
  }

  # Facet
  if(!is.null(group_by)) {
    if(!is(df[,group_by], "factor") & !is(df[,group_by], "character")){
      cat(paste0(ep, "group_by must be a categorical variable (factor or character). Returning un-facetted plot.\n"))
    } else {
          p = p + facet_wrap(vars(!!sym(group_by)))
    }
  }

  return(p)
}

#' @importFrom ggplot2 .data aes position_dodge ggplot geom_violin geom_boxplot guides
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line xlab scale_fill_manual
#' @importFrom qualpalr qualpal
#' @importFrom ggbeeswarm geom_quasirandom

.makeViolin <- function(df, x, y, plot_cells = FALSE, color_by = NULL, color_palette = NULL) {

  ep = "{cellula::.makeViolin() via plot_Coldata()} - "

  # Define mappings
  aes_cd = aes(y = .data[[y]])

  if(!is.null(x)) {
    if(!is(df[,x], "factor") & !is(df[,x], "character") & !is(df[,x], "logical"))
      stop(paste0(ep, "x must be a categorical variable (coercible to factor)"))
    df[,x] = factor(df[,x])
    aes_cd$x = aes(x = .data[[x]])$x
    if(is.null(color_by)) {
      aes_cd$colour = aes(colour = .data[[x]])$colour
      aes_cd$fill = aes(fill = .data[[x]])$fill
    }
    if(!is.null(color_by)) {
      if(!is(df[,color_by], "factor") & !is(df[,color_by], "character") & !is(df[,color_by], "logical"))
        stop(paste0(ep, "color_by must be a categorical variable (coercible to factor)"))
      aes_cd$colour = aes(colour = .data[[color_by]])$colour
      aes_cd$fill = aes(fill = .data[[color_by]])$fill
    }
  } else {
    df$group = y
    aes_cd$x = aes(x = .data[["group"]])$x
    if(!is.null(color_by)) {
      df[,color_by] = factor(df[,color_by])
      aes_cd$x = aes(x = .data[[color_by]])$x
      aes_cd$colour = aes(colour = .data[[color_by]])$colour
      aes_cd$fill = aes(fill = .data[[color_by]])$fill
    }
  }

  dodge <- position_dodge(width = 1)

  # Violins
  p <- ggplot(df, mapping = aes_cd) +
    geom_violin(alpha = 0.35,
                scale = "width",
                width = 0.9,
                position = dodge)

  # Beeswarm
  if(plot_cells) {
    if(is.null(x) & !is.null(color_by)) {
      aes_bee = aes_cd
      aes_bee$x = aes(x = .data[[color_by]])$x
      aes_bee$colour = aes(colour = .data[[color_by]])$colour
      pwidth = 1/length(levels(df[,color_by]))
    } else if(!is.null(x) & is.null(color_by)){
      pwidth = 1/length(levels(df[,x]))
      aes_bee = aes_cd
      aes_bee$x = aes(x = .data[[x]])$x
      aes_bee$colour = aes(colour = .data[[x]])$colour
    } else if(!is.null(x) & !is.null(color_by)){
      if(x != color_by) {
        aes_bee = aes(y = .data[[y]], x = .data[[x]])
        pwscale = length(unique(paste0(df[,x], "_", df[,color_by])))
        pwidth = 1/pwscale
      }

    }

    p <- p + geom_quasirandom(mapping = aes_bee,
                              groupOnX = TRUE,
                              #width = pwidth,
                              bandwidth = 1,
                              alpha = 0.3,
                              size = 0.6,
                              dodge.width = 1)
  }

  # Boxplot
  p <- p + geom_boxplot(mapping = aes_cd,
                        width = 0.1,
                        outlier.shape = NA,
                        col = "black",
                        position = dodge,
                        show.legend = FALSE) +
    theme_classic() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
    ) +
    xlab(NULL)

  # Colors
  if(!is.null(x) | !is.null(color_by)){
    if(!is.null(x) & is.null(color_by)) group_var = x else group_var = color_by
        if(is.null(color_palette)) {
          if(length(unique(df[,group_var])) > 1) {
            pal = qualpal(n = length(unique(df[,group_var])))$hex
          } else {
            pal = "gray"
          }

    } else {
      pal = color_palette
    }

    cscale = scale_colour_manual(values = pal)
    fscale = scale_fill_manual(values = pal)
    fguides = guides(boxplot = "none")
    p <- p + cscale + fscale + fguides
  }

  return(p)
}


#' @importFrom ggplot2 .data aes ggplot stat_density_2d after_stat
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line
#' @importFrom qualpalr qualpal
#' @importFrom colorspace scale_fill_continuous_sequential

.makeScatter <- function(df, x, y,
                         color_by = NULL,
                         color_palette = NULL,
                         contour = TRUE) {
  # Define mappings
  aes_cd = aes(x = .data[[x]], y = .data[[y]])

   if(!is.null(color_by)) {
      aes_cd$colour = aes(colour = .data[[color_by]])$colour
    }
  p <- ggplot(df, mapping = aes_cd) +
    geom_point(alpha = 1,
               size = 0.5)
  if(contour) {
    p <- p + stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(.data[["level"]])),
                    alpha = 0.5,
                    bins = 10,
                    col = "black",
                    linetype = 4,
                    linewidth = 0.2) +
      scale_fill_continuous_sequential(palette = "Heat 2")
  }
    p <- p + theme_classic() +
             theme(plot.margin = margin(1, 1, 1, 1, "cm"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               text = element_text(family = "sans")
        )

    # Colors
    if(!is.null(color_by)) {
      if(is.null(color_palette)) {
        pal = qualpal(n = length(unique(df[,color_by])))$hex
      } else {
        pal = color_palette
      }
      cscale = scale_colour_manual(values = pal)
      p = p + cscale
    }

  return(p)
}


#' @importFrom ggplot2 .data aes ggplot geom_tile element_text scale_x_discrete
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line scale_fill_gradientn
#' @importFrom stats dist hclust
#' @importFrom colorspace sequential_hcl

.makeHeatmap <- function(df, x, y, clustered = TRUE, color_palette = NULL) {

  # Calculate Jaccard index
  jdf = expand.grid(levels(df[,x]), levels(df[,y]))

  jdf$intersection = apply(jdf[,1:2], 1, function(row) {
    length(intersect(which(df[,1] == row[1]), which(df[,2] == row[2])))
  })

  jdf$union = apply(jdf[,1:2], 1, function(row) {
    length(union(which(df[,1] == row[1]), which(df[,2] == row[2])))
  })

  jdf$`Jaccard index`= jdf$intersection/jdf$union

  colnames(jdf)[c(1,2)] = c(x, y)

  # Reorder (protect against 1, 10, 2, ...)
  jdf[,x] = .reorderNumericLevels(jdf[,x])
  jdf[,y] = .reorderNumericLevels(jdf[,y], rev = TRUE)

  # Clustering
  if(clustered) {
    mat = matrix(jdf$`Jaccard index`, nrow = length(unique(jdf[,x])))
    col_ord = hclust(dist(mat))$order
    row_ord = rev(hclust(dist(t(mat)))$order)
    jdf[,x] = factor(jdf[,x], levels = unique(jdf[,x])[col_ord])
    jdf[,y] = factor(jdf[,y], levels = unique(jdf[,y])[row_ord])
  }

  # Plot
  p = ggplot(jdf, aes(x = .data[[x]], y = .data[[y]])) +
           geom_tile(aes(fill = .data[["Jaccard index"]])) +
           scale_x_discrete(position = "top")

  p = p + theme_classic() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          text = element_text(family = "sans")
    )

  # Colors
  if(is.null(color_palette)) {
    pal = rev(sequential_hcl(n = 40, palette = "YlGnBu"))[10:40]
  } else {
    pal = color_palette
  }

  cscale = scale_fill_gradientn(colours = pal)

  p = p + cscale
  return(p)

}


.reorderNumericLevels <- function(f, rev = FALSE) {
  conv = suppressWarnings(as.numeric(levels(f)))
    if(!any(is.na(conv))) {
    fl = as.character(sort(as.numeric(levels(f))))
    if(rev) fl = fl[rev(seq_along(fl))]
    levels(f) = fl
    return(f)
  } else {
    return(f)
  }
}

