#' Plot cluster modularity
#'
#' Plots a heatmap showing pairwise cluster modularity 
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param name character, the name in the \code{metadata slot} of the
#'     \code{SingleCellExperiment} object, e.g. "modularity_SNN_100"
#' @param type character, one of \code{"heatmap"} or \code{"graph"}
#'     for different plot types
#'
#' @return a heatmap of pairwise modularity
#'
#' @importFrom igraph graph_from_adjacency_matrix E layout_with_lgl
#' @importFrom ggplot2 scale_x_discrete geom_tile element_blank element_text
#' @importFrom ggplot2 ggplot .data theme_minimal
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export

plotModularity <- function(sce, name, type = "heatmap") {
    name <- paste0("modularity_", name)
    modmat <- log2(metadata(sce)[[name]] + 1)

    if (type == "heatmap") {
        m <- as.data.frame(expand.grid(
            seq_along(rownames(modmat)),
            seq_along(rownames(modmat))
        ))
        colnames(m) <- c("x", "y")
        m$value <- vapply(
            seq_len(nrow(m)),
            function(x) modmat[m[x, 1], m[x, 2]], numeric(1)
        )
        p <- ggplot(m, aes(x = .data[["x"]], y = .data[["y"]])) +
            geom_tile(aes(fill = .data[["value"]], color = NA)) +
            scale_x_discrete(position = "top")
        p <- p + theme_minimal() +
            theme(
                plot.margin = margin(1, 1, 1, 1, "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                text = element_text(family = "sans"),
                title = "log2(modularity ratio + 1)"
            )
    } else if (type == "graph") {
        cgr <- graph_from_adjacency_matrix(modmat,
            mode = "upper",
            weighted = TRUE,
            diag = FALSE
        )
        plot(cgr,
            edge.width = E(cgr)$weight * 5,
            layout = igraph::layout_with_lgl,
            main = name
        )
    }
}

#' Plot approximate silhouette widths
#'
#' Plots a beesswarm plot showing the approximate width of each cell within
#' each cluster
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param name character, the name in the metadata slot of the \code{SingleCellExperiment} object, e.g.
#'     "silhouette_SNN_100"
#'
#' @return a beeswarm plot of silhouette widths
#'
#' @importFrom ggplot2 ggplot aes_string theme_minimal
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom S4Vectors metadata metadata<-
#'
#' @details This plotting function has been directly taken from the OSCA book
#'     ("Orchestrating Single Cell Analysis with Bioconductor").
#'
#' @export

plotSilhouette <- function(sce, name) {
    name <- paste0("silhouette_", name)
    ggplot(metadata(sce)[[name]], aes_string(x = "cluster", y = "width", colour = "closest")) +
        geom_quasirandom(method = "smiley") +
        theme_minimal()
}

#' Plot proportions
#'
#' Plots the proportion of cells in each label and grouping
#'
#' @param sce a SingleCellExperiment object
#' @param label character, the column name in the
#'    colData of the SingleCellExperiment used for
#'    colouring the plot
#' @param group_by character, the column name in the
#'    colData of the SingleCellExperiment used for
#'    grouping the bars
#' @param color_palette named character, the name of the color
#'    palette to be used. Default is \code{NULL},
#'    meaning an automatic palette will be generated.
#'
#' @return a stacked bar plot of proportions of cells in each
#'    label, grouped according to the \code{group_by} variable.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_bw
#' @importFrom ggplot2 scale_y_continuous theme element_text element_blank
#' @importFrom ggplot2 geom_hline .data ylab
#' @importFrom grDevices colors
#'
#' @export

plotProportions <- function(sce,
                            label,
                            group_by,
                            color_palette = NULL) {
    ## Sanity checks
    # Error prefix
    ep <- .redm("{cellula::plotProportions()} - ")

    # Checks
    if (!is(sce, "SingleCellExperiment")) {
        stop(paste0(ep, "Must provide a SingleCellExperiment object"))
    }
    if (!label %in% colnames(colData(sce))) {
        stop(paste0(ep, "label was not found in the colData of the object"))
    }
    if (!group_by %in% colnames(colData(sce))) {
        stop(paste0(ep, "group_by was not found in the colData of the object"))
    }
    if (!is.null(color_palette)) {
        if (length(names(color_palette)) == 0) {
            stop(paste0(ep, "color_palette must be a named vector"))
        }
        if (!all(unique(colData(sce)[, label]) %in% names(color_palette))) {
            stop(paste0(ep, "color_palette does not contain names for all the levels of label"))
        }
    }

    if (is.null(color_palette)) {
        total_cols <- length(unique(colData(sce)[, label]))
        if (total_cols <= 30) {
            color_palette <- .cpal_qual(n = total_cols)
        } else {
            color_palette <- colors()[sample(seq_len(657), size = total_cols)]
        }
        names(color_palette) <- unique(as.character(colData(sce)[, label]))
    }

    cd <- as.data.frame(colData(sce)[, c(label, group_by), drop = FALSE])

    tab <- as.data.frame(table(cd))
    tab <- do.call(rbind, split(tab, tab[, group_by]))
    tab$prop <- unlist(lapply(split(tab, tab[, group_by]), function(x) x$Freq / sum(x$Freq)))

    ggplot(tab, aes(x = .data[[group_by]], y = .data[["prop"]] * 100, fill = .data[[label]])) +
        geom_bar(stat = "identity", color = "black") +
        geom_hline(yintercept = 0, linewidth = 0.5) +
        scale_fill_manual(values = color_palette) +
        theme_bw() +
        ylab("Proportion") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line.y = element_line()
        )
}

#' Plot Dimensionality Reduction
#'
#' Plots dimensionality reduction coordinates from a SingleCellExperiment object,
#' coloured and grouped according to colData fields
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param dr character, the name in the reducedDim slot of \code{sce}, e.g.
#'     \code{"UMAP"}
#' @param color_by character, column name in the \code{colData} slot of \code{sce},
#'     e.g. "cluster". Will be used to assign colors. Automatically detects
#'     whether the variable is categorical or continuous.
#' @param shape_by character, column name in the \code{colData} slot of  \code{sce},
#'     e.g. "cluster". Will be used to assign shapes. Can only be used for
#'     categorical variables.
#' @param group_by character, column name in the \code{colData} slot of the \code{sce},
#'     e.g. "individual". Will be used to facet the plot. Can only be used for
#'     categorical variables.
#' @param label_by character, column name in the \code{colData} slot of \code{sce},
#'     e.g. "cluster". Will be used to add labels to the plot. Can only be used
#'     for categorical variables.
#' @param knn_smooth logical, should the continuous values be smoothed using
#' 	   k-nearest neighbors? Default is FALSE.
#' @param smoothing_k numeric, the number of neighbors to use for smoothing.
#'     Default is 10.
#' @param point_size numeric, the size of the points in the plot. Default is 0.7.
#' @param label_size numeric, the size of the font for the labels. Default is 2.
#' @param plotting_order character, either "max", "min" or "random". This will determine
#'     the plotting order of the points in the plot. Default is "max", meaning that
#'     points will be plotted showing the highest values of \code{color_by} on top.
#'     If "min", points will be plotted in the order of the minimum value of \code{color_by}.
#'     If "random", points will be plotted in a random order.
#' @param outline logical, should a black outline be painted around the point cloud?
#'     Default is FALSE
#' @param outline_size numeric, the thickness of the outline, expressed as a fraction
#'     of the dot size. Default is 1.3, meaning the outline will be point size * 1.3.
#' @param arrows logical, should two perpendicular arrows be drawn on the bottom left
#'     corner of the plot be drawn? Default is TRUE.
#' @param exprs_use character, the name of the assay to be used when plotting a
#'     feature. Default is \code{"logcounts"}.
#' @param num_scale numeric or character, either "auto" (default), which computes
#'     color scale numbers automatically, or a numeric vector of length 2 with
#'     lower and upper limits for the scale.
#' @param color_palette a character string containing colors to be used. Default
#'     is \code{NULL}, meaning an automatic palette will be generated based on
#'     the type of datum supplied in \code{color_by}. Some palettes can be named:
#'     "Sunset", "Parula", "Turbo", "YlGnBU" for quantitative data, and "Qualpal",
#'     "Tritan", "Protan", "Tableau", "Pear", "Polychrome" and "Polylight" for
#'     categorical data. The "Protan" and "Tritan" palettes are adapted for CVD
#'     (Color Vision Deficiency) caused by protanopia and tritanopia respectively.
#' @param trajectory a character string indicating the `metadata` slot containing
#'     segment trajectories to be plotted. Usually either
#'     \code{"Slingshot_embedded_curves"} or \code{"Monocle_embedded_curves"}.
#' @param rescale logical, should coordinates be rescaled between 0 and 1? Default is TRUE.
#'
#' @return a ggplot object showing the dimensionality reduction coordinates colored, shaped and/or
#'     faceted as set by the arguments.
#'
#' @details This plotting function is heavily inspired by the \code{\link[scater]{plotReducedDim}}
#' function from the \code{{scater}} package.
#'
#' 	   The function can be used to plot 2D dimensionality reduction coordinates, coloured by
#'     a continuous or categorical variable, shaped by a categorical variable, and facetted
#'     by a categorical variable. It can also add labels to the plot, and draw arrows in the
#'     bottom left corner with coordinate names.
#'
#'     The function also allows to plot a feature/gene from the \code{assay} slot of
#'     the \code{SingleCellExperiment} object by looking up the character in the columns
#' 	   of the \code{rowData} slot.
#'
#'     To avoid overplotting the \code{knn_smooth} argument can be set to TRUE,
#'     which will smooth numeric values using k-nearest neighbors taking, for
#'     each cell, the mean of the values of its neighbourhood. The size of the
#' 	   neighborhood can be set with the \code{smoothing_k} argument.
#'
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes scale_colour_gradientn guides geom_point
#' @importFrom ggplot2 guide_legend .data theme_void geom_segment geom_text coord_fixed
#' @importFrom ggplot2 facet_wrap vars arrow unit geom_label scale_colour_manual theme
#' @importFrom ggplot2 margin theme_minimal
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats median complete.cases
#' @importFrom methods is
#'
#' @export

plot_DR <- function(sce,
                    dr = "UMAP",
                    color_by = NULL,
                    shape_by = NULL,
                    group_by = NULL,
                    label_by = NULL,
                    knn_smooth = FALSE,
                    smoothing_k = 10,
                    point_size = 0.7,
                    label_size = 2,
                    plotting_order = "max",
                    outline = FALSE,
                    outline_size = 1.3,
                    arrows = TRUE,
                    exprs_use = "logcounts",
                    num_scale = "auto",
                    color_palette = NULL,
                    trajectory = NULL,
                    rescale = TRUE) {
    ## Sanity checks
    # Error prefix
    ep <- .redm("{cellula::plot_DR()} - ")

    # Checks
    if (!is(sce, "SingleCellExperiment")) {
        stop(paste0(ep, "Must provide a SingleCellExperiment object"))
    }
    if (!dr %in% reducedDimNames(sce)) {
        stop(paste0(ep, "Could not find a slot named \"", dr, "\" in the SingleCellExperiment object"))
    }
    if (!is.null(color_by)) {
        if (!color_by %in% colnames(colData(sce)) & !color_by %in% rowData(sce)$Symbol & !color_by %in% rowData(sce)$ID) {
            stop(paste0(ep, color_by, " is not a column or a rowData element of the SingleCellExperiment object"))
        }
    }
    if (!is.null(shape_by)) {
        if (!shape_by %in% colnames(colData(sce))) {
            stop(paste0(ep, shape_by, " column not found in the colData of the SingleCellExperiment object"))
        }
    }
    if (!is.null(label_by)) {
        if (!label_by %in% colnames(colData(sce))) {
            stop(paste0(label_by, " column not found in the colData of the SingleCellExperiment object"))
        }
        if (!is(colData(sce)[, label_by], "character") & !is(colData(sce)[, label_by], "factor")) {
            stop(paste0(ep, "Cannot label by a numeric value, convert it to factor first."))
        }
    }
    if (!is.null(group_by)) {
        if (!group_by %in% colnames(colData(sce))) {
            stop(paste0(ep, group_by, " column not found in the colData of the SingleCellExperiment object"))
        }
        if (!is(colData(sce)[, group_by], "character") & !is(colData(sce)[, group_by], "factor")) {
            stop(paste0(ep, "Cannot group by a numeric value, convert it to factor first."))
        }
    }
    if (num_scale != "auto" & !is(num_scale, "numeric")) {
        stop(paste0(ep, " `num_scale` must be either a numeric vector or \"auto\""))
    }

    if (is(num_scale, "numeric") & length(num_scale) != 2) {
        stop(paste0(ep, " if specified, `num_scale` must contain exactly 2 numbers"))
    }

    udf <- as.data.frame(reducedDim(sce, dr)[, c(1, 2)])
    udf <- udf[complete.cases(udf), ]

    # Rescaling
    if (rescale) {
        orig_range <- lapply(udf, range)
        names(orig_range) <- c("x", "y")
        udf[, c(1, 2)] <- apply(udf[, c(1, 2)], 2, function(x) (x - min(x)) / diff(range(x)))
    }

    colnames(udf) <- c("x", "y")

    if (!is.null(color_by) & (color_by %in% rowData(sce)$Symbol | color_by %in% rowData(sce)$ID | color_by %in% rownames(sce)) &
        !color_by %in% colnames(colData(sce))) {
        feature <- which(rowData(sce)$Symbol == color_by | rowData(sce)$ID == color_by | rownames(sce) == color_by)
        colData(sce)[, color_by] <- as.numeric(assay(sce[feature, ], exprs_use))
    }

    # Check which mappings to include
    include <- list(color_by, shape_by, group_by, label_by)
    include <- unique(unlist(include[!is.null(include)]))

    for (i in include) udf[[i]] <- colData(sce)[rownames(udf), i]

    colnames(udf)[seq_len(length(include)) + 2] <- include

    classes <- unlist(lapply(udf[, include, drop = FALSE], class))

    names(classes) <- include
    aes_umap <- aes(x = .data[["x"]], y = .data[["y"]])

    # Add color aesthetic if needed and decide what type of scale
    if (!is.null(color_by)) {
        aes_umap$colour <- aes(colour = .data[[color_by]])$colour

        if (classes[color_by] == "numeric") {
            if (knn_smooth) {
                udf[, color_by] <- .smoothValues(udf[, c("x", "y")], k = smoothing_k)
            }

            if (plotting_order == "max") {
                udf <- udf[order(udf[, color_by]), ]
            } else if (plotting_order == "min") {
                udf <- udf[order(-udf[, color_by]), ]
            } else if (plotting_order == "random") {
                udf <- udf[sample(seq_len(nrow(udf)), nrow(udf)), ]
            }

            udf <- rbind(
                udf[is.na(udf[, color_by]), ],
                udf[!is.na(udf[, color_by]), ]
            )

            pal <- .choosePalette(cpal = color_palette, default = "Sunset")

            if (num_scale != "auto") {
                lims <- num_scale
                cscale <- scale_colour_gradientn(
                    colours = pal, na.value = "lightgray",
                    limits = lims
                )
            } else {
                cscale <- scale_colour_gradientn(colours = pal, na.value = "lightgray")
            }
            cguides <- NULL
        } else if (classes[color_by] %in% c("character", "logical")) {
            udf[, color_by] <- factor(udf[, color_by])
            if (length(levels(udf[, color_by])) == 1) {
                pal <- .choosePalette(cpal = NA)
            } else {
                pal <- .choosePalette(
                    cpal = color_palette, default = "qualpal",
                    n = length(unique(udf[, color_by]))
                )
            }
            if (!is.null(names(pal)) & all(names(pal) %in% levels(udf[, color_by]))) {
                pal <- pal[levels(udf[, color_by])]
            }
            cscale <- scale_colour_manual(values = pal, na.value = "lightgray")
            cguides <- guides(color = guide_legend(override.aes = list(size = 2)))
        } else if (classes[color_by] == "factor") {
            pal <- .choosePalette(
                cpal = color_palette, default = "qualpal",
                n = length(unique(udf[, color_by]))
            )

            if (!is.null(names(pal)) & all(names(pal) %in% levels(udf[, color_by]))) {
                pal <- pal[levels(udf[, color_by])]
            }
            cscale <- scale_colour_manual(values = pal, na.value = "lightgray")
            cguides <- guides(color = guide_legend(override.aes = list(size = 2)))
        }
    }

    # Add shape aesthetic if needed
    if (!is.null(shape_by)) {
        if (classes[shape_by] == "numeric") {
            stop("Cannot map `shape` aesthetic to a numeric value")
        } else if (classes[shape_by] == "character") {
            udf[, shape_by] <- factor(udf[, shape_by])
        }
        aes_umap$shape <- aes(shape = .data[[shape_by]])$shape
    }

    if (!is.null(label_by)) {
        labels <- unique(udf[, label_by])
        medoids <- lapply(labels, function(x) {
            t(apply(udf[udf[, label_by] == x, c("x", "y")], 2, median))
        })
        label_df <- as.data.frame(do.call(rbind, medoids))
        label_df[, label_by] <- labels
        colnames(label_df) <- c("x", "y", "label")
    }

    if (!is.null(group_by)) {
        if (is(udf[, group_by], "character")) {
            udf[, group_by] <- factor(udf[, group_by])
        }
    }

    # Plot construction
    p <- ggplot(udf, mapping = aes_umap)

    if (outline) {
        p <- p + geom_point(
            size = point_size * outline_size,
            color = "black"
        )
    }

    p <- p + geom_point(size = point_size) +
        theme(plot.margin = margin(2, 2, 2, 2, "cm"))

    # Coordinate arrows (bottom left corner)

    if (arrows) {
        arrow_1 <- data.frame(x = 0, xend = 0.15, y = 0, yend = 0)
        arrow_2 <- data.frame(x = 0, xend = 0, y = 0, yend = 0.15)
        text_1 <- data.frame(x = -0.03, y = 0.02)
        text_2 <- data.frame(x = 0.02, y = -0.03)

        p <- p +
            geom_segment(
                data = arrow_1,
                mapping = aes(
                    x = .data[["x"]],
                    y = .data[["y"]],
                    xend = .data[["xend"]],
                    yend = .data[["yend"]]
                ),
                inherit.aes = FALSE,
                color = "black",
                linewidth = 0.3,
                arrow = arrow(
                    length = unit(0.2, "cm"),
                    type = "closed"
                )
            ) +
            geom_segment(
                data = arrow_2,
                mapping = aes(
                    x = .data[["x"]],
                    y = .data[["y"]],
                    xend = .data[["xend"]],
                    yend = .data[["yend"]]
                ),
                color = "black",
                inherit.aes = FALSE,
                linewidth = 0.3,
                arrow = arrow(
                    length = unit(0.2, "cm"),
                    type = "closed"
                )
            ) +
            geom_text(
                data = text_1,
                aes(
                    x = .data[["x"]],
                    y = .data[["y"]]
                ),
                inherit.aes = FALSE,
                label = paste0(dr, " 2"),
                angle = 90,
                size = 2,
                hjust = "left",
                vjust = "top",
                color = "black"
            ) +
            geom_text(
                data = text_2,
                aes(
                    x = .data[["x"]],
                    y = .data[["y"]]
                ),
                inherit.aes = FALSE,
                label = paste0(dr, " 1"),
                size = 2,
                hjust = "left",
                vjust = "bottom",
                color = "black"
            )
    }

    p <- p + coord_fixed()

    if (!is.null(color_by)) {
        p <- p + cscale + cguides
    }

    if (!is.null(label_by)) {
        p <- p + geom_label_repel(
            data = label_df,
            mapping = aes(
                x = .data[["x"]],
                y = .data[["y"]],
                label = .data[["label"]]
            ),
            inherit.aes = FALSE, size = label_size
        )
    }

    if (!is.null(group_by)) {
        p <- p + facet_wrap(group_by)
    }

    if (!is.null(trajectory)) {
        trj <- metadata(sce)[[trajectory]]

        if (rescale) {
            trj$x0 <- .rescalen(trj$x0,
                from = range(orig_range$x),
                to = range(udf$x)
            )
            trj$x1 <- .rescalen(trj$x1,
                from = range(orig_range$x),
                to = range(udf$x)
            )

            trj$y0 <- .rescalen(trj$y0,
                from = range(orig_range$y),
                to = range(udf$y)
            )

            trj$y1 <- .rescalen(trj$y1,
                from = range(orig_range$y),
                to = range(udf$y)
            )
        }
        trj <- as.data.frame(trj)
        p <- p + geom_segment(
            data = trj,
            mapping = aes(
                x = .data[["x0"]],
                xend = .data[["x1"]],
                y = .data[["y0"]],
                yend = .data[["y1"]]
            ),
            inherit.aes = FALSE
        )
    }
    p +
        theme_minimal() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
        )
}

#' Plot a gene expression dot plot
#'
#' Plots a dot plot with gene expression from a SingleCellExperiment object
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param genes character string, the genes to be plotted (matched with rownames)
#' @param group_by character, column name in the \code{colData} slot of \code{sce},
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
  ep <- .redm("{cellula::plot_dots()} - ")

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
  if(!exprs_use %in% assayNames(sce))
    stop(paste0(ep, "Could not find an assay named \"", exprs_use, "\" in the SingleCellExperiment object"))
  if(!is(colData(sce)[,group_by], "character") & !is(colData(sce)[,group_by], "factor"))
      stop(paste0(ep, "Cannot group by a numeric value, convert it to factor first."))
  if(any(!genes %in% rownames(sce))){
    found <- intersect(genes, rownames(sce))
    if(length(found) == 0) {
      stop(paste0(ep, "None of the genes were found in the SingleCellExperiment object"))
    } else {
      message(paste0(ep, "Warning: some genes (", (length(genes) -length(found)), ") were not found in the dataset."))
      genes <- found
    }
  }

  sce_byclust <- lapply(unique(colData(sce)[,group_by]), function(x) {
    cur <- assay(sce, exprs_use)[genes, which(colData(sce)[,group_by] == x)]
    props <- apply(cur, 1, function(x) sum(x > threshold))/ncol(cur)
    aves <- apply(cur, 1, function(x) mean(x))
    final_df <- data.frame("proportion" = props,
                          "mean_expression" = aves,
                          "cluster" = x,
                          "gene" = genes)
    return(final_df)
  })

  if(cluster_genes) {
    scdf_props_genes <- do.call(cbind, lapply(sce_byclust, function(x) x$proportion))
    scdf_exp_genes <- do.call(cbind, lapply(sce_byclust, function(x) x$mean_expression))
    hc_props_genes <- hclust(dist(scdf_props_genes))$order
    hc_exp_genes <- hclust(dist(scdf_exp_genes))$order
  }
  
  if(cluster_groups) {
    scdf_props_clusters <- do.call(rbind, lapply(sce_byclust, function(x) x$proportion))
    scdf_exp_clusters <- do.call(rbind, lapply(sce_byclust, function(x) x$mean_expression))
    hc_props_clusters <- hclust(dist(scdf_props_clusters))$order
    hc_exp_clusters <- hclust(dist(scdf_exp_clusters))$order
  }
  
  scdf <- do.call(rbind, sce_byclust)
  
  if(cluster_genes) {
    scdf$gene <- factor(scdf$gene, levels = genes[as.numeric(hc_exp_genes)])
  } else scdf$gene <- factor(scdf$gene, levels = rev(genes))
  if(cluster_groups) {
    scdf$cluster <- factor(scdf$cluster, levels = unique(colData(sce)[,group_by])[as.numeric(hc_exp_clusters)])
  } 

  scdf$mean_expression[scdf$mean_expression == 0] <- NA
  scdf$proportion[scdf$proportion == 0] <- NA

  if(is.null(color_palette)) {
    pal <- .cpal_seq_ylgnbu()
  } else {
    pal <- color_palette
  }
  cscale <- scale_colour_gradientn(colours = pal, na.value = "lightgray")

  if(format == "wide") {
    p <- ggplot(scdf, aes(y = .data[["cluster"]],
                         x = .data[["gene"]],
                         size = .data[["proportion"]],
                         colour = .data[["mean_expression"]])
               ) +
      geom_point(shape = 16, na.rm = TRUE) +
      cscale +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
            ) +
      labs(colour = "Mean expression", size = "Proportion")
  } else if(format == "tall") {
    p <- ggplot(scdf, aes(y = .data[["gene"]],
                         x = .data[["cluster"]],
                         size = .data[["proportion"]],
                         colour = .data[["mean_expression"]])
               ) +
      geom_point(shape = 16, na.rm = TRUE) +
      cscale +
      theme_minimal() +
      theme(axis.text.y = element_text(face = "italic"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
            ) +
      labs(colour = "Mean expression", size = "Proportion")
  }
  p
}

#' @importFrom scater retrieveCellInfo

.makeColdataDF <- function(sce, columns) {
  do.call(cbind, lapply(columns, function(x) {
    ret <- retrieveCellInfo(sce, by = x)
    col <- data.frame(ret$value, row.names = colnames(sce))
    colnames(col) <- ret$name
    return(col)
  })
  )
}

#' Plot cell metadata
#'
#' Plot data from the colData of the SCE object
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param y character, the column from \code{colData(sce)} whose values need to be plotted
#' @param x character, the column from \code{colData(sce)} to be used as a grouping variable
#'     default is \code{NULL}, which means all points belong to one violin plot only
#' @param color_by character, the column name from \code{colData(sce)} to be used as a
#'     colouring variable. If \code{NULL} (default), violins (and points) will not be colored.
#' @param group_by character, the column name from \code{colData(sce)} to be used as a
#'     facetting variable. Default is \code{NULL}.
#' @param color_palette character string of colors for categorical plots or heatmap.
#'     Default is \code{NULL}.
#' @param contour logical, should contours be plotted on the scatterplot?
#'     Default is \code{TRUE}.
#' @param clustered logical, should rows and columns be clustered in the confusion
#'     matrix heatmap? Default is \code{TRUE}.
#' @param plot_cells logical, should cells be plotted as well? Default is \code{FALSE}
#'     for speed.
#' @param rotate_labels logical, should x axis labels be rotated by 45 degrees 
#'     to prevent overlapping text? Only applies to heatmaps and violin plots. 
#'     Default is \code{TRUE}     
#'
#' @returns a ggplot with different types of visualization depending on the
#'     classes of the columns chosen. 
#'     
#' @details If y is a numeric and x is NULL, or a categorical (character/factor), 
#'     the function plots a violin + boxplot divided by x and/or by color_by. 
#'     If y and x are both numeric, the function plots a scatterplot with optional 
#'     2D kernel contouring. If y and x are both categorical, the function plots 
#'     a heatmap of the confusion matrix, showing pairwise Jaccard indices. 
#'     All plots can be facetted.
#'
#' @importFrom ggplot2 facet_wrap
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
                         plot_cells = FALSE,
                         rotate_labels = TRUE) {
  ## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::plot_Coldata()} - ")
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))

  # Prepare data
  include <- c(x, y, color_by, group_by)
  include <- unique(include[!is.null(include)])

  if(any(!include %in% colnames(colData(sce))))
    stop(paste0(ep, "Some columns were not found in the colData"))

  df <- .makeColdataDF(sce, include)

  if(!is.null(x) & is(df[,x], "numeric") & is(df[,y], "numeric")){
    p <- .makeScatter(df, x, y, color_by, color_palette, contour)
  } else if((!is(df[,x], "numeric") & !is.null(x)
             & is(df[,y], "numeric"))
            | (is.null(x) & is(df[,y], "numeric"))) {
    p <- .makeViolin(df, x, y, plot_cells, color_by, color_palette, rotate_labels)
  } else if(!is.null(x)
            & (is(df[,x], "character") | is(df[,x], "factor"))
            & (is(df[,y], "character") | is(df[,y], "factor"))) {
    df[,x] <- factor(df[,x])
    df[,y] <- factor(df[,y])
    p <- .makeHeatmap(df, x, y, clustered, color_palette, rotate_labels)
  }

  # Facet
  if(!is.null(group_by)) {
    if(!is(df[,group_by], "factor") & !is(df[,group_by], "character")){
      cat(paste0(ep, "group_by must be a categorical variable (factor or character). Returning un-facetted plot.\n"))
    } else {
          p <- p + facet_wrap(group_by)
    }
  }
  p
}


#' Plot a stacked violin plot
#'
#' Plots a stacked violin plot of at least two features in different groupings
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param features character vector, list of features from \code{sce}
#' @param cluster character, column name in the \code{colData} slot of \code{sce},
#'     e.g. "cluster". Will be used to assign colors. 
#'     Must be a categorical variable.
#' @param split_by character, column name in the \code{colData} slot of the \code{sce},
#'     e.g. "individual". Will be used to divide violin plots. 
#'     Must be a categorical variable.
#' @param stack logical, should all violin plots be stacked on 1 column? 
#'     Default is TRUE.
#' @param de (optional) data frame containing the results of differential expression.
#'     Standardized column names are enforced.
#' @param significance logical, should a bar with adjusted p-value and log2(FC) values
#'     be plotted on top of the violin plots? Default is FALSE. Requires \code{de} to
#'     be specified as well.
#' @param color_palette character, a vector of colors to be used for the \code{fill}
#'     aesthetic for the violin plot areas.
#' @param nrow numeric, the number of rows to build the facetted plot. Only used
#'     if \code{stacked = FALSE}.
#' @param ncol numeric, the number of columns to build the facetted plot. Only used
#'     if \code{stacked = FALSE}.
#'
#' @return a ggplot object showing a facetted violin plot with several features grouped
#'     and colored according to two categorical variables.
#' 
#' @details This function has two possible outcomes, controlled by the \code{stack} argument: 
#'     1) \code{stack = TRUE} returns a stacked violin plot where each feature is a row, and 
#'     violins correspond to whatever grouping is specified in "cluster"; 
#'     2) \code{stack = FALSE} returns a faceted grid of violin plots where each feature is a
#'     panel.
#' 
#'     The first behaviour is useful to check markers/gene expression patterns in different 
#'     groupings with an idea of the distribution. It reproduces \code{scanpy}'s stacked
#'     violin plot. 
#'     The second behaviour is useful to show comparisons between groupings/conditions 
#'     for a set of features, and allows differential expression significance information
#'     to be overlaid. 
#' 
#'     Since the idea is to make different groupings and features comparable
#'     the \code{\link[ggplot2]{geom_violin}} geom is set with \code{scale = "width"}, which can show the
#'     distribution of the data but does not scale with the number of data points (cells).
#'
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom ggplot2 ggplot aes .data theme_minimal theme element_text element_blank
#' @importFrom ggplot2 geom_violin geom_segment geom_text element_line
#' @importFrom ggplot2 facet_wrap vars scale_fill_manual ggplot_build
#' @importFrom reshape2 melt
#' @importFrom methods is
#'
#' @export

stacked_Violin <- function(sce, features, cluster, split_by = NULL, stack = TRUE,
                          de = NULL, significance = FALSE, color_palette = NULL,
                          nrow = NULL, ncol = NULL) {

## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::stacked_Violin()} - ")

  # Checks
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(is.null(features))
      stop(paste0(ep, "Must provide features to plot"))
  if(length(features) < 2)
      stop(paste0(ep, "Must provide at least two features to plot. To plot a single feature use `plot_Coldata()`"))     
  if(is.null(cluster))
      stop(paste0(ep, "Must provide a cluster column name"))
  if(!cluster %in% colnames(colData(sce)))
      stop(paste0(ep, cluster, " is not a column of the SingleCellExperiment object"))
  if(!is(colData(sce)[,cluster], "character") & !is(colData(sce)[,cluster], "factor"))
      stop(paste0(ep, "`cluster` cannot be a numeric value, convert it to factor first."))

  if(!is.null(split_by)) {
    if(!split_by %in% colnames(colData(sce)))
      stop(paste0(ep, split_by, " is not a column of the SingleCellExperiment object"))
    if(!is(colData(sce)[,split_by], "character") & !is(colData(sce)[,split_by], "factor"))
      stop(paste0(ep, " `split_by` cannot be a numeric value, convert it to factor first."))
  }
  if(!is.null(color_palette)){
    if(length(color_palette) < length(unique(colData(sce)[,cluster])))
      stop(paste0(ep, "the color palette cannot be shorter than the number of unique values of `cluster`")) 
  }

  if(significance & is.null(de))
      stop(paste0(ep, "`signifiance` cannot be TRUE if you do not provide `de`."))

  if(is.null(split_by)) split_by = cluster
  
    feat_index <- which(rowData(sce)$Symbol %in% features | rowData(sce)$ID %in% features | rownames(sce) %in% features)
   
    featuredf = as.data.frame(as.matrix(t(assay(sce, "logcounts")[feat_index,])))
    colnames(featuredf) = features = rowData(sce)[colnames(featuredf), "Symbol"]
    
    featuredf$cluster = colData(sce)[,cluster]
    featuredf$splitby = colData(sce)[,split_by]
    fdf_melt = reshape2::melt(featuredf, value.name = "expression")
    fdf_melt$variable = factor(fdf_melt$variable, levels = colnames(featuredf))
    
    if(stack) {
      p = ggplot(fdf_melt, aes(y = .data[["expression"]], 
                               x = .data[["cluster"]], 
                               fill = .data[["splitby"]])) + 
          geom_violin(scale = "width") + 
        facet_wrap(~factor(variable, levels = features),
                  nrow = length(features), 
                  strip.position = "right", 
                  scales = "free",
                  axis.labels = "all_y") +
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_line(),
              #axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15)) 

      if(!is.null(color_palette)) {
        p = p + scale_fill_manual(values = color_palette)
    }
      
    } else {
      p = ggplot(fdf_melt, aes(y = .data[["expression"]], 
                               x = .data[["cluster"]], 
                               fill = .data[["splitby"]])) + 
        geom_violin(scale = "width") + 
        geom_boxplot(data = fdf_melt, mapping =  aes(y = .data[["expression"]], x = .data[["cluster"]]),
                     width = 0.2,
                     position = position_dodge(width = 0.9), outliers = FALSE) +
        ylim(0, max(fdf_melt$expression) * 1.6) +
        facet_wrap(~factor(variable, levels = features), nrow = nrow, ncol = ncol) +
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.line = element_line(),
              axis.ticks = element_line(),
              #axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12)) 
      if(!is.null(color_palette)) {
        p = p + scale_fill_manual(values = color_palette)
      }
    }
    
    if(significance) {
      
      built = ggplot_build(p)$data[[1]]
      xcoords = unique(built$x)
      ymax = max(built$y)*1.1
      
      sigdf = as.data.frame(de[de$gene %in% features, c("gene", "logFC", "padj")])
      sigdf$variable = sigdf$gene
      sigdf$logFC = round(sigdf$logFC, digits = 2)
      sigdf$padj = formatC(sigdf$padj, format = "e", digits = 2)
      sigdf$padj[sigdf$padj == "0.00e+00"] = "< 1e-299"
      sigdf$lab = paste0("logFC: ", sigdf$logFC, "\n", "p: ", sigdf$padj)
      
      p = p + geom_segment(aes(x = xcoords[1], xend = xcoords[2], 
                               y = ymax, yend = ymax), 
                            linewidth = 0.5) +
          geom_text(data = sigdf, 
                   aes(x = min(xcoords) + 0.45, y = ymax * 1.3, label = .data[["lab"]]), 
                   inherit.aes = FALSE, nudge_x = -.225, size = 2.5) 
      
    }
    p
}


#' Plot a multi-panel dimensionality reduction 
#'
#' Plots a multi-panel dimensionality reduction where points are colored 
#' by different features (one per panel)
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param dr character, the name of the dimensional reduction slot (retrieved through 
#'     \code{reducedDim(sce, dr)}). Default is "UMAP".
#' @param dims numeric, vector of 2 dimensions to plot. Default is 1, 2
#' @param features character vector with features (e.g. genes) from \code{sce}. 
#' @param point_size numeric, the size of the points in the plot. Default is 1.2
#' @param plot_order character, one of "decreasing" (default), "increasing", or "random". 
#'    Influences the way points are plotted, important when dealing with overplotting.
#' @param exprs_use character, the name of the \code{assay} in the object whose values will
#'    be plotted. Default is "logcounts".
#' @param rng_seed numeric, the random number generator seed used when \code{plot_order = "random"}.
#' @param common_scale logical, should the points have a single, common color scale (TRUE) or 
#'    should each panel have its own scale? Deafult is FALSE. 
#' @param color_palette character, vector of colors to be interpolated across for the color
#'    aesthetic. Default is NULL meaning a standard quantitative palette from \code{cellula}
#'    will be used.
#' @return a ggplot object showing a multi-panel plot of different features - one feature per panel -
#'    where there is either a common scale or an individual scale per panel.
#' 
#' @details This function is inspired by Seurat's \code{\link[Seurat]{FeaturePlot}}. 
#'     This function has two possible outcomes, controlled by the \code{common_scale} argument: 
#'     1) \code{common_scale = TRUE} returns a faceted plot where there is a single scale for color.
#'     This can be useful to compare gene expression across genes using color. 
#'     2) \code{common_scale = FALSE} returns a \code{patchwork} grid of feature plots, each with its own
#'     individual scale. 
#' 
#' @importFrom SummarizedExperiment rowData assay 
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot aes .data theme_minimal theme element_blank
#' @importFrom ggplot2 geom_point facet_wrap 
#' @importFrom ggplot2 scale_colour_gradientn coord_fixed
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @importFrom patchwork wrap_plots
#'
#' @export


multipanel_DR <- function(sce, dr = "UMAP", dims = c(1,2), 
                          features, point_size = 1.2, 
                          plot_order = "decreasing",
                          exprs_use = "logcounts",
                          rng_seed = 11,
                          common_scale = FALSE,
                          color_palette = NULL) {

 ## Sanity checks
  # Error prefix
  ep <- .redm("{cellula::multipanel_DR()} - ")

  # Checks                          
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(!dr %in% reducedDimNames(sce))
    stop(paste0(ep, "the `dr` provided was not found in this object."))
  if(is.null(features))
      stop(paste0(ep, "Must provide features to plot"))
  if(length(features) < 2)
      stop(paste0(ep, "Must provide at least two features to plot. To plot a single feature use `plot_DR()`"))     
  if(!plot_order %in% c("decreasing", "increasing", "random"))
      stop(paste0(ep, "`plot_order` must be one of \"decreasing\", \"increasing\", or \"random\""))
  if(!is(rng_seed, "numeric"))
      stop(paste0(ep, "`rng_seed` must be a numeric"))

  feat_index <- which(rowData(sce)$Symbol %in% features | rowData(sce)$ID %in% features | rownames(sce) %in% features)
  
  featuredf = as.data.frame(as(t(assay(sce, exprs_use)[feat_index,]), "matrix"))
  colnames(featuredf) = features = rowData(sce)[colnames(featuredf), "Symbol"]
  
  fdf_melt = reshape2::melt(featuredf, value.name = "expression")
  fdf_melt$variable = factor(fdf_melt$variable, levels = colnames(featuredf))
  fdf_melt$x = rep(reducedDim(sce, dr)[,dims[1]], length(features))
  fdf_melt$y = rep(reducedDim(sce, dr)[,dims[2]], length(features))
  
  if(is.null(color_palette)) {
    colorpal = .cpal_seq_ylgnbu()
  } else {
    colorpal = color_palette
  }
  
  if(plot_order == "decreasing") {
    fdf_melt = fdf_melt[order(fdf_melt$expression),]
  } else if(plot_order == "increasing") {
    fdf_melt = fdf_melt[order(fdf_melt$expression, decreasing = TRUE),]
  } else if(plot_order == "random") {
    set.seed(rng_seed)
    new_ord = sample(seq_len(nrow(fdf_melt)), size = nrow(fdf_melt), replace = FALSE)
    fdf_melt = fdf_melt[new_ord,]
  }
  
  if(common_scale) {
    p = ggplot(fdf_melt, aes(x = .data[["x"]], 
                             y = .data[["y"]], 
                             color = .data[["expression"]])) + 
      geom_point(size = point_size) + 
      facet_wrap(~factor(variable, levels = features)) +
      theme_minimal() + 
      theme(panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()) +
      scale_colour_gradientn(colours = colorpal) +
      coord_fixed()
  } else {
    plotlist = lapply(split(fdf_melt, fdf_melt$variable), function(x) {
      ggplot(x, aes(x = .data[["x"]], 
                    y = .data[["y"]], 
                    color = .data[["expression"]])) + 
        geom_point(size = point_size) + 
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank()) +
        scale_colour_gradientn(colours = colorpal, name = unique(x$variable)) +
        coord_fixed()
    })
      p = patchwork::wrap_plots(plotlist = plotlist)
  }
  p
}

#' @importFrom ggplot2 .data aes position_dodge ggplot geom_violin geom_boxplot guides
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line xlab scale_fill_manual
#' @importFrom ggplot2 element_text
#' @importFrom ggbeeswarm geom_quasirandom

.makeViolin <- function(df, x, y, plot_cells = FALSE, color_by = NULL, color_palette = NULL,
                        rotate_labels = TRUE) {

  ep <- .redm("{cellula::.makeViolin() via plot_Coldata()} - ")

  # Define mappings
  aes_cd <- aes(y = .data[[y]])

  if(!is.null(x)) {
    if(!is(df[,x], "factor") & !is(df[,x], "character") & !is(df[,x], "logical"))
      stop(paste0(ep, "x must be a categorical variable (coercible to factor)"))
    df[,x] <- factor(df[,x])
    aes_cd$x <- aes(x = .data[[x]])$x
    if(is.null(color_by)) {
      aes_cd$colour <- aes(colour = .data[[x]])$colour
      aes_cd$fill <- aes(fill = .data[[x]])$fill
    }
    if(!is.null(color_by)) {
      if(!is(df[,color_by], "factor") & !is(df[,color_by], "character") & !is(df[,color_by], "logical"))
        stop(paste0(ep, "color_by must be a categorical variable (coercible to factor)"))
      aes_cd$colour <- aes(colour = .data[[color_by]])$colour
      aes_cd$fill <- aes(fill = .data[[color_by]])$fill
    }
  } else {
    df$group <- y
    aes_cd$x <- aes(x = .data[["group"]])$x
    if(!is.null(color_by)) {
      df[,color_by] <- factor(df[,color_by])
      aes_cd$x <- aes(x = .data[[color_by]])$x
      aes_cd$colour <- aes(colour = .data[[color_by]])$colour
      aes_cd$fill <- aes(fill = .data[[color_by]])$fill
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
      aes_bee <- aes_cd
      aes_bee$x <- aes(x = .data[[color_by]])$x
      aes_bee$colour <- aes(colour = .data[[color_by]])$colour
      #pwidth <- 1/length(levels(df[,color_by]))
    } else if(!is.null(x) & is.null(color_by)){
      #pwidth <- 1/length(levels(df[,x]))
      aes_bee <- aes_cd
      aes_bee$x <- aes(x = .data[[x]])$x
      aes_bee$colour <- aes(colour = .data[[x]])$colour
    } else if(!is.null(x) & !is.null(color_by)){
      if(x != color_by) {
        aes_bee <- aes(y = .data[[y]], x = .data[[x]])
        pwscale <- length(unique(paste0(df[,x], "_", df[,color_by])))
        #pwidth <- 1/pwscale
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
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          plot.margin = margin(1, 1, 1, 1, "cm"),
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
            pal <- .cpal_qual(n = length(unique(df[,group_var])))
          } else {
            pal <- "gray"
          }

    } else {
      pal <- color_palette
    }

    cscale <- scale_colour_manual(values = pal, na.value = "lightgray")
    fscale <- scale_fill_manual(values = pal, na.value = "lightgray")
    fguides <- guides(boxplot = "none")
    p <- p + cscale + fscale + fguides 
  }
  
  if(rotate_labels) {
    p <- p + theme(axis.text.x.bottom = element_text(angle = 45, 
                                                     hjust = 1))
  }
  p
}

#' @importFrom ggplot2 .data aes ggplot stat_density_2d after_stat scale_fill_gradientn
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line

.makeScatter <- function(df, x, y,
                         color_by = NULL,
                         color_palette = NULL,
                         contour = TRUE) {
  # Define mappings
  aes_cd <- aes(x = .data[[x]], y = .data[[y]])

   if(!is.null(color_by)) {
      aes_cd$colour <- aes(colour = .data[[color_by]])$colour
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
      scale_fill_gradientn(colours = .cpal_seq_heat2())
  }
    p <- p + theme_minimal() +
             theme(plot.margin = margin(1, 1, 1, 1, "cm"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"),
               text = element_text(family = "sans")
        )

    # Colors
    if(!is.null(color_by)) {
      if(is.null(color_palette)) {
        pal <- .cpal_qual(n = length(unique(df[,color_by])))
      } else {
        pal <- color_palette
      }
      cscale <- scale_colour_manual(values = pal, na.value = "lightgray")
      p <- p + cscale
    }

  return(p)
}


#' @importFrom ggplot2 .data aes ggplot geom_tile element_text scale_x_discrete
#' @importFrom ggplot2 theme_minimal theme element_blank labs element_line scale_fill_gradientn
#' @importFrom stats dist hclust

.makeHeatmap <- function(df, x, y, clustered = TRUE, color_palette = NULL, 
                         rotate_labels = TRUE) {

  if(!is(df[,x], "factor")) df[,x] = factor(df[,x])
  if(!is(df[,y], "factor")) df[,y] = factor(df[,y])
  
  # Calculate Jaccard index
  jdf <- expand.grid(levels(df[,x]), levels(df[,y]))

  jdf$intersection <- apply(jdf[,c(1,2)], 1, function(row) {
    length(intersect(which(df[,1] == row[1]), which(df[,2] == row[2])))
  })

  jdf$union <- apply(jdf[,c(1,2)], 1, function(row) {
    length(union(which(df[,1] == row[1]), which(df[,2] == row[2])))
  })

  jdf$`Jaccard index`= jdf$intersection/jdf$union

  colnames(jdf)[c(1,2)] <- c(x, y)

  # Reorder (protect against 1, 10, 2, ...)
  #jdf[,x] <- .reorderNumericLevels(jdf[,x])
  #jdf[,y] <- .reorderNumericLevels(jdf[,y], rev = TRUE)

  # Clustering
  if(clustered) {
    mat <- matrix(jdf$`Jaccard index`, nrow = length(unique(jdf[,x])))
    col_ord <- hclust(dist(mat))$order
    row_ord <- rev(hclust(dist(t(mat)))$order)
    jdf[,x] <- factor(jdf[,x], levels = unique(jdf[,x])[col_ord])
    jdf[,y] <- factor(jdf[,y], levels = unique(jdf[,y])[row_ord])
  }

  # Plot
  p <- ggplot(jdf, aes(x = .data[[x]], y = .data[[y]])) +
           geom_tile(aes(fill = .data[["Jaccard index"]])) +
           scale_x_discrete(position = "top")

  p <- p + theme_minimal() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          text = element_text(family = "sans")
    )

  if(rotate_labels) p <- p + theme(axis.text.x.top = element_text(angle = 45, 
                                                                 hjust = 0))
  # Colors
  if(is.null(color_palette)) {
    pal <- .cpal_seq_ylgnbu()
  } else {
    pal <- color_palette
  }
  cscale <- scale_fill_gradientn(colours = pal)
  p = p + cscale
  p
}

#' Gene expression heatmap
#' 
#' Draws a heatmap of gene expression from a SingleCellExperiment object
#' 
#' @param sce a SingleCellExperiment object
#' @param genes character, the genes to include in the heatmap.
#'     These will be looked up in the rownames of the object, which are also the 
#' 	   rownames of the data frames in the DE results. At least 2 genes must be provided.
#' @param coldata_cols character, the columns in the \code{colData} slot of the \code{sce} object
#'     to include in the heatmap as annotations. Their order influences column ordering
#' 	   when \code{cluster_cols = FALSE}. Default is \code{NULL}.
#' @param order_by character, the columns in the \code{colData} slot of the \code{sce} object
#'     to order the heatmap by. Default is \code{coldata_cols}. Overrides \code{cluster_cols} if
#'     not \code{NULL}.
#' @param scale logical, should the expression data be scaled by row? Default is \code{TRUE}
#' @param exprs character, the expression values to use. Default is \code{"logcounts"}
#' @param aggregate logical, should the data be aggregated? If \code{TRUE}, the data will be
#'     aggregated by the variables indicdated in \code{aggregate_by}. Default is \code{FALSE}
#' @param aggregate_by character, the columns in the \code{colData} slot of the \code{sce} object
#'     to aggregate by. Default is \code{NULL}
#' @param aggregate_fun character, the function to use for aggregation. Default is \code{"sum"}
#' @param color_pal character, the color palette to use for the heatmap. Default is \code{NULL}
#' @param annotation_pal named list of named character vectors, the color palette 
#' 	   to use for the annotations. Default is \code{NULL}. See Details
#' @param cluster_cells logical, should the columns (cells) be clustered? 
#'     Default is \code{TRUE}. It is automatically set to \code{FALSE} if \code{order_by} is 
#'     not \code{NULL}.
#' @param cluster_genes logical, should the rows (genes) be clustered? Default is \code{TRUE}
#' @param gaps character, the columns in \code{coldata_cols} to leave gaps between. See Details
#'     Default is \code{NULL}.
#' @param clip_values numeric, the percentile to clip the values at. Default is \code{99.5}
#' @param raster logical, should the heatmap be rasterized? Default is \code{FALSE}, but 
#'     it is changed to \code{TRUE} automatically if the number of cells or genes is 
#' 	   greater than 2000. This behaviour can be overriden by setting \code{raster = "force"}.
#' @param ... additional arguments to pass to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return a heatmap of the gene expression data, optionally with annotations 
#'    from the \code{colData} slot
#' 
#' @details this function generates a heatmap of genes using expression data from the \code{sce} 
#' 	  object (specified by \code{exprs}) with optional aggregation. The heatmap can be annotated
#' 	  with any number of column annotations, i.e. colouring segments above each column according
#'    to the \code{colData} slot of the \code{sce} object. 
#' 		
#' 	  The \code{annotation_pal} argument is a named list of named character (colour) vectors, 
#' 	  where each element of the list corresponds to a column in \code{coldata_cols}. 
#' 	  The named character vectors should have the same length as the number of unique values 
#'    in the column, and the names should be the unique values. 
#' 	
#' 	  The order of colData columns in \code{order_by} influences the order of the columns in the
#' 	  heatmap. If \code{cluster_cols = FALSE}, the columns will be ordered according to the order
#' 	  of the columns in \code{order_by}, where the first element dictates the major order, then
#' 	  the second element the minor order, etc, e.g. if \code{order_by = c("condition", "label")}
#'    columns will be first ordered by condition and then, within each condition, by label.
#' 	  If \code{cluster_cols = TRUE}, the columns will be ordered according to the clustering 
#' 	  of the columns. Two numeric columns cannot be properly used for ordering since they will
#' 	  have an individual level for each cell.
#' 
#' 	  The aggregation arguments \{code{aggregate}, \code{aggregate_by}, \code{aggregate_fun} 
#'    are used to control if and how the data is aggregated. If \code{aggregate = TRUE}, the data
#'    will be aggregated by the variables in \code{aggregate_by} using the function specified by 
#'    \code{aggregate_fun}. 
#' 
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom stats quantile
#' @importFrom grDevices colorRampPalette colors
#' 
#' @export
 
plotHeatmap <- function(sce, 
                        genes, 
                        coldata_cols = NULL, 
						order_by = coldata_cols,
                        scale = TRUE,
						exprs = "logcounts",
                        aggregate = FALSE,
                        aggregate_by = NULL,
						aggregate_fun = "sum",
                        color_pal = NULL, 
                        annotation_pal = NULL,
                        cluster_cells = TRUE,
                        cluster_genes = TRUE,
                        gaps = NULL,
                        clip_values = 99.5,
                        raster = FALSE,
						...){
  
  ### Sanity checks
	# error prefix
	ep = .redm("cellula::plotHeatmap() - ")
	
  dependencies = data.frame("package" = c("ComplexHeatmap", "circlize"),
                            "repo" = c("BioC", "CRAN"))
  if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(is.null(genes))
    stop(paste0(ep, "Must provide genes to plot"))
  if(length(genes) < 2)
    stop(paste0(ep, "Must provide at least two genes to plot"))
  if(!is.null(coldata_cols)) {
    if(!all(coldata_cols %in% colnames(colData(sce))))
      stop(paste0(ep, "Some columns were not found in the colData"))
	if(!is.null(order_by)) {
		if(!all(order_by %in% colnames(colData(sce))))
			stop(paste0(ep, "Some columns for ordering were not found in the colData"))
		if(!all(order_by %in% coldata_cols))
			stop(paste0(ep, "Some columns for ordering were not found in coldata_cols"))	
	}
    if(!is.null(gaps)) {
      if(!gaps %in% coldata_cols)
        stop(paste0(ep, "the gaps column must be included in coldata_cols"))
    }
	if(!is(clip_values, "numeric"))
		stop(paste0(ep, "clip_values must be a numeric value"))
	if(clip_values <= 0 | clip_values > 100)
		stop(paste0(ep, "clip_values must be between 0 excluded and 100 included"))	
  }
  
  feat_index <- which(rowData(sce)$Symbol %in% genes | rowData(sce)$ID %in% genes | rownames(sce) %in% genes)
  
  if(length(feat_index) == 0)
    stop(paste0(ep, "No genes were found in the object"))
  
  if(aggregate) {
    if(is.null(aggregate_by))
      stop(paste0(ep, "aggregate_by must not be NULL if aggregate is TRUE"))
    if(!all(aggregate_by %in% colnames(colData(sce))))
      stop(paste0(ep, "Some columns for aggregation were not found in the colData"))
    
    agg_vector = apply(colData(sce)[,aggregate_by,drop=FALSE], 1, function(x) paste(x, collapse = "__"))
    agg = aggregateAcrossCells(sce,
                               ids = agg_vector,
                               statistics = aggregate_fun, 
                               use.assay.type = "counts")
		if(exprs == "logcounts") {
			agg = computeLibraryFactors(agg)
			agg = logNormCounts(agg)
		}
		mat = assay(agg, exprs)[feat_index,]
		cd = as.data.frame(colData(agg))
		non_na_cols = colnames(cd)[unlist(lapply(cd, function(x) !all(is.na(x))))]
		order_by = intersect(order_by, non_na_cols)
		coldata_cols = intersect(coldata_cols, non_na_cols)
	} else {
		mat = as(assay(sce, exprs)[feat_index,], "matrix") #oh god forgive me
		cd = as.data.frame(colData(sce))
	}

	mat[is.na(mat)] = 0
	remove = which(rowSums(mat) == 0)

	if(length(remove) > 0) {
		mat = mat[-remove,]
		message("Removed ", length(remove), " genes with no expression")
	}

  if(scale) mat = t(scale(t(mat)))

  if(!is.null(coldata_cols) & !is.null(order_by)) {
	if(length(order_by) > 1) {
    	ordered_cols = do.call(order, cd[,order_by])
	} else {
		ordered_cols = order(cd[,order_by])
	}
    cd = cd[ordered_cols,]
    mat = mat[,ordered_cols]

    if(is.null(annotation_pal)){	

    	categorical_cols_test = unlist(lapply(cd[,coldata_cols,drop=FALSE], function(x) 
      		inherits(x, "character")|inherits(x, "factor")|inherits(x,"logical")))      
    	categorical_cols = categorical_cols_test[categorical_cols_test]
		
		numeric_cols_test = unlist(lapply(cd[,coldata_cols,drop=FALSE], function(x) 
      		inherits(x, "numeric")|inherits(x, "integer")))
		numeric_cols = numeric_cols_test[numeric_cols_test]

		# Palette for categorical columns
		if(length(categorical_cols) > 0) {

			num_cat_cols = lapply(names(categorical_cols), function(x) {
								length(unique(as.character(cd[,x,drop=TRUE])))
								})
			
			names(num_cat_cols) = names(categorical_cols)
			total_cols = sum(unlist(num_cat_cols))
			
			if(total_cols <= 30) {
				pal_auto = .cpal_qual(n = total_cols) 
			} else {
				pal_auto = colors()[sample(seq_len(657), size = total_cols)]
			}
			if(length(num_cat_cols) > 1) {
				indices = c(0, cumsum(unlist(num_cat_cols)))
				col_list_categorical = lapply(seq_len(length(indices)-1), function(i) {
						pal = pal_auto[(indices[i]+1):(indices[i+1])]
						names(pal) = unique(as.character(cd[,names(categorical_cols)[i]]))
						pal
					})
			} else {
				names(pal_auto) = unique(as.character(cd[,names(categorical_cols)]))
				col_list_categorical <- list(pal_auto)
				
			}
				names(col_list_categorical) = names(categorical_cols)
		} else {
			col_list_categorical = NULL
		}
		
		# Palette for numeric columns
		if(length(numeric_cols) > 0) {
			cols_choose = .cpal_qual_pear(length(numeric_cols))
			col_list_numeric = lapply(seq_len(length(numeric_cols)), function(i) {
				circlize::colorRamp2(breaks = c(0, max(cd[,names(numeric_cols)[i]], na.rm = TRUE)), 
				colors = c("white", cols_choose[i]))
			})
			names(col_list_numeric) = names(numeric_cols)
		} else {
			col_list_numeric = NULL
		}
		col_list = c(col_list_categorical, col_list_numeric)
		col_list = col_list[!unlist(lapply(col_list, is.null))]
  	} else {
    	col_list = annotation_pal
  	}
	# Final annotation
	column_ha = ComplexHeatmap::HeatmapAnnotation(df = cd[,names(col_list), drop=FALSE],
													  col = col_list)
  }  
  if(!raster | raster != "force")
    if(nrow(mat) > 20000 | ncol(mat) > 20000) raster = TRUE
  if(raster == "force") raster = FALSE

  if(!is.null(gaps)) {
    gaps = cd[,gaps]
  }
	# Color mapping and clipping for the heatmap
  if(scale) {
    if(!is.null(clip_values)) {
      clip_values = clip_values*10
      cr = max(abs(quantile(mat, seq(0,1,length.out=1001), na.rm = TRUE)[c(1001-clip_values, clip_values)]))
      
    } else {
      cr = max(abs(range(mat[!is.na(mat)])))
    }
    color_range = c(-cr, 0, cr)
  } else {
    if(!is.null(clip_values)){
	# Clip to percentiles with variable decimal digit precision
	  decimal_places = nchar(strsplit(as.character(clip_values), "\\.")[[1]][2])
      if(is.na(decimal_places)) decimal_places = 0
	  clip_values = clip_values * 10^(decimal_places)
      cr = quantile(mat, seq(0, 1, length.out = (10^(decimal_places+2)+1)), na.rm = TRUE)[clip_values]
      color_range = range(mat[!is.na(mat)])
      color_range[which.max(color_range)] = cr
    } else {
      color_range = range(mat[!is.na(mat)])
    }
  }
  
  if(!is.null(color_pal)) {
      color_pal_ramp = colorRampPalette(color_pal)(25)
  } else {
    if(scale) {
      color_pal_ramp = colorRampPalette(c("blue", "white", "red"))(25)
    } else {
      color_pal_ramp = colorRampPalette(.cpal_seq_ylgnbu())(25)
    }
  }
  
  col_fun = circlize::colorRamp2(quantile(color_range, 
                                          probs = seq(0, 1, length.out = 25)), 
                                 		  color_pal_ramp)
  
  transf_name = ifelse(scale, "Z score", exprs)

  if(aggregate) {
	name_scale = paste0(transf_name, "\n(", exprs, "(", aggregate_fun, "))")
  } else {
	name_scale = transf_name
  }

  if(!is.null(coldata_cols)){
    H = ComplexHeatmap::Heatmap(mat, 
                                row_names_gp = grid::gpar(cex = 1), 
                                column_names_gp = grid::gpar(cex = 1), 
                                col = col_fun,
                                top_annotation = column_ha,
                                cluster_columns = cluster_cells,
                                cluster_rows = cluster_genes,
                                column_split = gaps,
                                use_raster = raster,
                                name = name_scale,
                                show_column_names = FALSE,
								...)
    
  } else {
    H = ComplexHeatmap::Heatmap(mat, 
                                row_names_gp = grid::gpar(cex = 1), 
                                column_names_gp = grid::gpar(cex = 1), 
                                col = col_fun,
                                cluster_columns = cluster_cells,
                                cluster_rows = cluster_genes,
                                use_raster = raster,
                                name = name_scale,
                                show_column_names = FALSE,
								...)
    
  }
  ComplexHeatmap::draw(H)
}

#' Pseudotime heatmap
#' 
#' Draws a heatmap of pseudotime-ordered gene expression 
#' 
#' @param sce a SingleCellExperiment object
#' @param genes character, the genes to include in the heatmap.
#' @param pseudotime character, the column in the \code{colData} 
#'     slot of the \code{sce} object with pseudotime values.
#' @param exprs character, the name of the assay with expression values to use. 
#'     Default is \code{"logcounts"}
#' @param cluster_genes logical, should the genes be clustered? 
#'     Default is \code{FALSE}
#' @param labels character, the column in the \code{colData} slot of the \code{sce} object
#'    to use as labels for heatmap annotation. Default is \code{NULL}
#' @param labels_pal named list of named character vectors, the color palette
#'    to use for the labels. Default is \code{NULL}
#' @param smoothing_window numeric, the window size for smoothing the expression values.
#'     Default is \code{200}
#'
#' @return a heatmap of the gene expression data ordered by pseudotime with,
#' 	   optionally, the labels of each cells as annotation on top of the heatmap.
#' 	   Each gene is scaled between 0 and 1 and smoothed using a zero-padded moving 
#' 	   average with the specified window size.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom stats median
#' 
#' @export 
plotPseudotimeHeatmap <- function(sce, 
							  genes, 
							  pseudotime, 
							  exprs = "logcounts", 
							  cluster_genes = FALSE,
                              labels = NULL, 
							  labels_pal = NULL, 
							  smoothing_window = 200) {

  ep = .redm("cellula::plotPseudotimeHeatmap() - ")
  dependencies = data.frame("package" = c("ComplexHeatmap", "circlize"),
                            "repo" = c("BioC", "CRAN"))
  if(checkFunctionDependencies(dependencies)) stop(paste0(ep, "Missing required packages."))
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(is.null(genes))
    stop(paste0(ep, "Must provide genes to plot"))
  if(length(genes) < 2)
    stop(paste0(ep, "Must provide at least two genes to plot"))
  if(!all(pseudotime %in% colnames(colData(sce))))
    stop(paste0(ep, "The pseudotime column was not found in the colData"))
  if(!all(labels %in% colnames(colData(sce))))
    stop(paste0(ep, "The labels column was not found in the colData"))
  if(!is(pseudotime, "character"))
    stop(paste0(ep, "The pseudotime column must be a character (column of colData)"))
  if(!is(smoothing_window, "numeric"))
    stop(paste0(ep, "The smoothing window must be a numeric"))
  if(!is.null(labels) & !is.null(labels_pal)) {
  	if(!all(names(labels_pal) %in% colData(sce)[,labels]))
   		stop(paste0(ep, "Not all levels in the labels_pal were found in the labels column"))
  }
  if(smoothing_window > length(colData(sce)[!is.na(colData(sce)[,pseudotime]),pseudotime]))
    stop(paste0(ep, "The smoothing window cannot be bigger than the number of non-NA pseudotime values"))

  cells_ok = which(!is.na(colData(sce)[,pseudotime]))
  sce = sce[,cells_ok]

  genemat = assay(sce, exprs)[genes, order(colData(sce)[,pseudotime])]
  
  genemat = t(apply(genemat, 1, function(x) 
    .rescalen(.mav(x, window = smoothing_window), 
                    to = c(0,1))))
  
  colors_pt = circlize::colorRamp2(colors = c("gray70", "purple"), 
                                   breaks = c(0, max(colData(sce)[,pseudotime])))
  
  if(!is.null(labels)) {
    labs = unique(colData(sce)[,labels])
    pt_means = lapply(labs, function(x) median(colData(sce)[colData(sce)[,labels] == x, pseudotime]))
    names(pt_means) = labs
    pt_means = sort(unlist(pt_means))
    
    if(is.null(labels_pal)) {
      labels_pal = .cpal_qual(length(labs))
      names(labels_pal) = sort(labs)
    }
    
    coldf = as.data.frame(colData(sce))
    coldf[,labels] = factor(coldf[,labels], 
                         levels = names(pt_means))
    
    lab_ord = colData(sce)[,labels][order(colData(sce)[,pseudotime])]
    lab_ord_mat = sapply(names(pt_means), function(x) as.numeric(lab_ord == x) + 1)
    
    df = as.data.frame(cbind(pseudotime = sort(colData(sce)[,pseudotime]), as.data.frame(lab_ord_mat)))
    
    col_list_group = lapply(names(pt_means), function(x) {
      cols = c("white", labels_pal[x])
      names(cols) = c("1", "2")
      cols
    })
    
    names(col_list_group) = names(pt_means)
    
    col_list = c(list(pseudotime = colors_pt), 
                 col_list_group) 
    
    column_ha = ComplexHeatmap::HeatmapAnnotation(df = df,
                                  col = col_list,
                                  show_legend = c(TRUE, rep(FALSE, length(unique(colData(sce)[,labels])))))
    
  } else {
    
    col_list = list(pseudotime = colors_pt)
    df = data.frame(pseudotime = sort(colData(sce)[,pseudotime]))
    column_ha = ComplexHeatmap::HeatmapAnnotation(df = df,
                                  				  col = col_list)
  }

  ComplexHeatmap::Heatmap(genemat,
						  col = .cpal_seq_parula(),
                          show_column_names = FALSE,
                          top_annotation = column_ha,
                          cluster_rows = cluster_genes,
                          cluster_columns = FALSE,
                          name = "Scaled expression")
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
			ep = .redm("cellula::plotDEHeatmap() - ")
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
	pal_cluster = .cpal_qual(length(dge))
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
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme_bw geom_hline geom_vline
#' @importFrom stats prcomp
#' 
#' @export 

plotLabelPCA <- function(sce, 
						replicates,
						labels,
						condition){

		### Sanity checks
		# error prefix	
			ep = .redm("cellula::plotLabelPCA() - ")

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

#' Per-label Mean-Difference plot
#' 
#' Plot a Mean-Difference plot of the DE results per label
#' 
#' @param dge a data frame with the differential expression results
#'    from the \code{doPBDGE} function
#' @param ntop numeric, the number of top genes to show in the plot
#' @param alpha numeric, the FDR threshold to use for filtering
#' @param lfc numeric, the log-fold change threshold to use for filtering
#' 
#' @return a faceted plot where each panel shows the mean (log(CPM)) vs logFC 
#'    for all genes in a label that do not have NA DE test results.
#' 
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap theme_bw geom_hline geom_vline
#' @importFrom ggrepel geom_text_repel
#' 
#' @export 

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

#' plot Receiver-Operator Characteristic curves
#' 
#' plot ROC curves with other classification statistics
#' 
#' @param sce a SingleCellExperiment object
#' @param gene character, the gene to use for the ROC curve
#' @param label character, the column in the \code{colData} slot of the \code{sce} object
#'   to use as labels for classification and facetting
#' @param exprs character, the name of the assay with expression values to use.
#'    Default is \code{"logcounts"}
#' @param plot logical, should the ROC curves be plotted? Default is \code{TRUE}
#' 
#' @return if \code{plot = TRUE}, a facetted plot with ROC curves 
#'    for each label in the \code{label} column together with other 
#'    classification statistics (see details). 
#' 	  If \code{plot = FALSE}, a list with the ROC curve data and other
#'    statistics is returned.
#' 
#' @details the function calculates the ROC curve for the gene specified 
#'     in the \code{gene} argument, showing the classification performance 
#' 	   of the gene for the label against all other labels. The function
#' 	   calculates the Area Under the Curve (AUC), Youden's J statistic,
#' 	   the entropy of the gene expression distribution, and the Jensen-Shannon
#' 	   Divergence between the gene expression distribution and an idealized 
#'     distribution where the gene is only expressed in one label.
#' 
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ggplot2 ggplot geom_path geom_abline facet_wrap theme_bw
#' @importFrom ggplot2 .data aes 
#' @importFrom stats embed
#' 
plotROC <- function(sce, 
					gene,
					label, 
					exprs = "logcounts", 
					plot = TRUE){

	rocs = lapply(unique(colData(sce)[,label]), function(x) {
		gin = colData(sce)[,label] == x
		truepos <- gin[order(assay(sce, exprs)[gene,], decreasing = FALSE)]
		TPR <- rev((sum(truepos) - cumsum(truepos))/sum(truepos))
		FPR <- rev(1-(cumsum(!truepos)/sum(!truepos)))

		# AUC
		h <- diff(FPR)
		ab = rowSums(embed(TPR, 2))
		auc = sum((ab*h)/2)
		message("AUC: ", round(auc, 4))
		
		# Youden's J: max(sens + spec - 1)
		# FPR = 1 - spec
		# spec = 1 - FPR
		# sens + spec - 1 = TPR + 1 - FPR - 1
		# sens + spec = TPR - FPR
		
		youden = which.max(TPR - FPR)
		thresh = sort(assay(sce, exprs)[gene, colData(sce)[,label] == x], decreasing = TRUE)
		youden_j = thresh[youden]
		
		# Entropy and Jensen-Shannon Divergence
		entropy = function(x) {
			x <- x[x > 0]
			x <- x/sum(x)
			-sum(x * log(x, 2))
		}
		
		real = unlist(lapply(split(assay(sce, "logcounts")[gene,], colData(sce)[,label]), mean))
		
		ideal = rep((1/(length(real)-1))/(length(real)-1), length(real))
		names(ideal) = names(real)
		ideal[x] = 1-1/(length(real)-1)
		ent = entropy(real/sum(real))
		mixture = ((real/sum(real)) + ideal) * 0.5
		jsd = entropy(mixture) - ((entropy(real) + entropy(ideal))*0.5)
		
		roc_curve = cbind(FPR, TPR)
		colnames(roc_curve) = c("FPR", "TPR")
		list(roc_curve = roc_curve,
			auc = auc,
			youden_j = youden_j,
			ent = ent,
			jsd = jsd)
	})	

		rocdf = as.data.frame(do.call(rbind, lapply(rocs, function(x) x$roc)))
		rocdf$label = rep(unique(colData(sce)[,label]), unlist(lapply(rocs, function(y) nrow(y$roc))))
		
		if(!plot) return(rocs)

	# Plot
	if(plot) { 
			ggplot(rocdf, aes(x = .data[["FPR"]], 
							  y = .data[["TPR"]], 
							  color = .data[["label"]])) +
				geom_path() +
				geom_abline(intercept = 0, slope = 1, linetype = 2) +
				facet_wrap("label") +
				theme_bw() 
	}
 
}	


#' Plot UMAP
#' 
#' Alias for plot_DR for consistency with previous iterations
#' 
#' @param ... see \code{\link{plot_DR}}
#' 
#' @return see \code{\link{plot_DR}}
#' 
#' @export 
#' 
 plot_UMAP <- function(...){
	plot_DR(...)
 }