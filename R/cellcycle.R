#' Assign cell cycle 
#' 
#' Calculate the cell cycle phase for every single cell using tricycle
#' 
#' @param sce a SingleCellExperiment object
#' @param recenter_method character, one of  \code{"optim"} (default),  \code{"chull"}, or  \code{"none"}
#' @param species character, one of  \code{"human"} (default) or  \code{"mouse"}.
#'    See \code{?tricycle::project_cycle_space()} for more information.
#' @param gname.type character, one of  \code{"SYMBOL"} (default) or \code{"ENSEMBL"}.
#'    Will be used to look up the genes in the  \code{rownames(sce)}.
#'    See \code{?tricycle::project_cycle_space()} for more information.
#' @param gname character, alternative \code{rowData(sce)} column with gene identifiers.
#'    Default is NULL. 
#'    see \code{?tricycle::project_cycle_space()} for more information.    
#' @param exprs_values character, name of the assay in \code{sce} to be used for
#'    projection. Default is \code{"logcounts"}.  
#' 
#' @returns a \code{SingleCellExperiment} object with the \code{tricycleEmbedding} 
#'    slot in \code{reducedDim(sce)} and \code{tricyclePosition} and 
#'    \code{tricyclePhase} columns in \code{colData(sce)}
#' 
#' @details Wrapper to the \code{tricycle} method (Zheng et al. 2022) for cell 
#'     cycle estimation, which uses a dataset with a continuum, periodic cell 
#'     cycle and projects the user data onto the embeddings. This wrapper only 
#'     adds two options for an automatic recentering procedure: \code{"optim"} 
#'     is an optimization routine to identify the point in the center which 
#'     is farthest from all other points, and \code{"chull"} calculates the center
#'     of the convex hull of the 2D embedding.
#'     
#'     Moreover, a crude phase assignment is calculated based on the angle (this
#'     is not part of the original publication but these assignments are detailed
#'     in the vignette.)
#'     
#' @importFrom SingleCellExperiment reducedDim 
#' @importFrom stats optim
#' @importFrom utils head installed.packages
#' @export

assignCellCycle <- function(sce, recenter_method = "optim", species = "human", 
                            gname.type = "SYMBOL", gname = NULL, 
                            exprs_values = "logcounts") {
  # Checks
  ep = .redm("{cellula::assignCellCycle()} - ")
  if (!"tricycle" %in% rownames(installed.packages()))
    stop(ep, "the `tricycle` package must be installed first.\n
                Run `BiocManager::install(\"tricycle\") to use this function.")
  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if (!species %in% c("human", "mouse")){
    stop(ep, species, " not recognized. Choose between \"human\" and \"mouse\"")
  }
  # Projection
  sce = tricycle::project_cycle_space(sce, gname.type = gname.type, gname = gname,
                                      species = species, exprs_values = exprs_values)
  embedding = reducedDim(sce, "tricycleEmbedding")
  # Recentering
  nc = switch(recenter_method, 
              "chull" = .findTricycleCenterChull(embedding),
              "optim" = .findTricycleCenterOpt(embedding),
              "none" = c(0,0))
  # Estimation of position
  sce = tricycle::estimate_cycle_position(sce, center.pc1 = nc[1], center.pc2 = nc[2])
  # Phase assignment according to the tricycle vignette
  phases = c(0, 0.25*pi, pi/2, pi, 1.5*pi, 1.75*pi, 2*pi)
  assigned = cut(sce$tricyclePosition, breaks = phases)
  levels(assigned) = c("G1/G0", "G1-S", "S", "G2-M", "M-G1", "G1/G0")
  sce$tricyclePhase = assigned
  sce
}


#' Plot cell cycle 
#' 
#' Plots the cell cycle assignment in a circular embedding
#' 
#' @param sce a SingleCellExperiment object
#' @param rings_by character, column in \code{colData(sce)} to divide the embedding by
#' @param color_by character, column in \code{colData(sce)}. Default is \code{"tricyclePosition"}.
#' @param cyclic_color character, one of "Cycle 1" or "Cycle 2" for a cyclic palette.
#'     Default is "Cycle 1"
#' @param color_palette character, any other color palette. Only used 
#' 
#' @returns a plot of cell cycle assignments
#' 
#' @importFrom SummarizedExperiment colData 
#' @importFrom ggplot2 ggplot aes .data geom_label geom_point scale_color_gradientn
#' @importFrom ggplot2 scale_colour_manual geom_path
#' @importFrom ggplot2 element_blank theme theme_minimal coord_fixed xlim ylim 
#' @importFrom methods is
#' 
#' @export

plotCycle <- function(sce, 
                      rings_by = NULL, 
                      color_by = "tricyclePosition",
                      cyclic_color = "Cycle 1",
                      color_palette = NULL) {
  
  # Checks
  ep = .redm("{cellula::plotCycle()} - ")
  if (!is(sce, "SingleCellExperiment"))
    stop(ep, "Must provide a SingleCellExperiment object")
  if(!is.null(color_by)){
    if (!color_by %in% colnames(colData(sce)))
      stop(ep, color_by, " is not a column of the SingleCellExperiment object")
  }
  if (!is.null(rings_by)){
    if (!rings_by %in% colnames(colData(sce)))
      stop(ep, rings_by, " is not a column of the SingleCellExperiment object")
  }
  if (!"tricyclePosition" %in% colnames(colData(sce)))
    stop(ep, "Tricycle positions not found. Run `assignCellCycle()` first.")
  
  # Colors
  cycol = match.arg(cyclic_color, c("Cycle 1", "Cycle 2"))
  # Create rings
  circlist = list()
  if(!is.null(rings_by)){
    rings = colData(sce)[,rings_by]
    nring = length(unique(rings))
  } else {
    rings = 1
    nring = 1
  }
  # Calculate positions on rings
  for(i in seq_along(unique(rings))) {
    cl = unique(rings)[i]
    curr = sce[,rings == cl]
    x = (cos(curr$tricyclePosition) * i * 0.5) 
    y = (sin(curr$tricyclePosition) * i * 0.5) 
    x[x>0] = x[x>0]+abs(rnorm(sum(x>0), 0, 0.04))
    x[x<0] = x[x<0]-abs(rnorm(sum(x<0), 0, 0.04))
    y[y>0] = y[y>0]+abs(rnorm(sum(y>0), 0, 0.04))
    y[y<0] = y[y<0]-abs(rnorm(sum(y<0), 0, 0.04))
    circlist[[i]] = data.frame(x = x, y = y, cluster = i)
    if(!is.null(color_by)) {
      circlist[[i]][,color_by] = colData(curr)[,color_by]
    } 
  }
  
  # Phase indicators
  phases = c(pi/2, pi, 1.5*pi, 0.25*pi, 1.75*pi)
  phases_x = cos(phases) * nring/2 * 1.3
  phases_y = sin(phases) * nring/2 * 1.3
  phase_labels = c(0.75*pi, 1.25*pi, 1.65*pi, 0, 0.35*pi)
  phase_labels_x = cos(phase_labels) * nring/2 * 1.4
  phase_labels_y = sin(phase_labels) * nring/2 * 1.4
  # DF for plotting
  dd = as.data.frame(do.call(rbind,circlist))
  if(is.null(color_by)) {
    color_by = "cells"
    dd[,"cells"] = "all"
  }
  # Circles
  circs = .makeCircles(r = seq(0.5, nring/2, by = 0.5))
  p = ggplot(dd, aes(x = .data[["x"]], y = .data[["y"]], col = .data[[color_by]]))
  for(i in unique(circs$r)){
    p = p + geom_path(data = circs[circs$r == i,], 
                      mapping = aes(x = x, y = y), 
                      inherit.aes = FALSE,
                      lwd = 0.15, 
                      lty = 2)
    }
    p = p + geom_segment(data = data.frame(rep(0,5), rep(0,5)), 
                         aes(x = 0, y = 0, xend = phases_x, yend = phases_y),
                         lwd = 0.2,
                         inherit.aes = FALSE)
  if(!is.null(rings_by)) {
    p = p + geom_label(data = data.frame(rep(0, nring), rep(0,nring)), 
                       mapping = aes(x = rep(0, nring), y = seq(0.5, nring/2, by = 0.5)-0.15, label = unique(rings)),
                       size = 2.5,
                       inherit.aes = FALSE)
  }
    p = p + 
    geom_text(data = data.frame(rep(0, 5), rep(0, 5)),
              aes(x = phase_labels_x, y = phase_labels_y, label = c("S", "G2-M", "M-G1", "G1/G0", "G1-S")),
              size = 5,
              inherit.aes = FALSE) +
    xlim(c(-1, 1) * nring/2 * 1.5) + 
    ylim(c(-1, 1) * nring/2 * 1.5) + 
    geom_point(alpha = 0.5)
    if(color_by == "tricyclePosition") {
      colorpal = .choosePalette(default = cycol)
      p = p + scale_color_gradientn(colors = colorpal)
    } else {
      if(is(dd[,color_by], "factor") | is(dd[,color_by], "character") | is(dd[,color_by], "logical")) {
        colorpal = .choosePalette(cpal = color_palette, n = length(unique(dd[,color_by])))
        p = p + scale_colour_manual(values = colorpal)
      } else {
        colorpal = .choosePalette(cpal = color_palette, default = "Sunset")
        p = p + scale_color_gradientn(colors = colorpal)
      }
    }
    p = p + coord_fixed() + 
            theme_minimal() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  axis.text.x = element_blank(), axis.text.y = element_blank(),
                  axis.title.x = element_blank(), axis.title.y = element_blank())
    p 
}

#' @noRd
.makeCircles <- function(r = 1, npt = 100){
  thetas <- seq(0, 2*pi, length.out = npt)
  circles = list()
  for(i in seq_along(r)) {
    xt <- r[i] * cos(thetas)
    yt <- r[i] * sin(thetas)
    circles[[i]] = data.frame(x = xt, y = yt, r = r[i])
  }
  as.data.frame(do.call(rbind, circles))
}

#' @noRd
#' @author Pierre-Luc Germain
.findTricycleCenterOpt <- function(e, nearest.q=.025, minImbalance=0.1){
  e <- as.matrix(t(e))
  nearest.q <- ceiling(nearest.q*length(ncol(e)))
  cf <- function(par){
    # ensure that we are not outside the cloud:
    if (!all(apply(e-par,1,FUN=function(x){
      min(sum(x>0),sum(x<=0))/length(x)>minImbalance
    }))) return(Inf)
    # get the sum of eucl dist to nearest cells
    d <- sqrt(colSums((e-par)^2))
    -sum(head(sort(d),nearest.q))
  }
  optim(c(x=0,y=0),cf)$par
}

#' @noRd
#' @importFrom grDevices chull
.findTricycleCenterChull <- function(e) {
  colMeans(e[chull(e),])
}