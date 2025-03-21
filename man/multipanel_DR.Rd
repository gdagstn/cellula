% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{multipanel_DR}
\alias{multipanel_DR}
\title{Plot a multi-panel dimensionality reduction}
\usage{
multipanel_DR(
  sce,
  dr = "UMAP",
  dims = c(1, 2),
  features,
  point_size = 1.2,
  plot_order = "decreasing",
  exprs_use = "logcounts",
  knn_smooth = FALSE,
  smoothing_k = 10,
  rng_seed = 11,
  common_scale = FALSE,
  color_palette = NULL
)
}
\arguments{
\item{sce}{a \code{SingleCellExperiment} object}

\item{dr}{character, the name of the dimensional reduction slot (retrieved through 
\code{reducedDim(sce, dr)}). Default is "UMAP".}

\item{dims}{numeric, vector of 2 dimensions to plot. Default is 1, 2}

\item{features}{character vector with features (e.g. genes) from \code{sce}.}

\item{point_size}{numeric, the size of the points in the plot. Default is 1.2}

\item{plot_order}{character, one of "decreasing" (default), "increasing", or "random". 
Influences the way points are plotted, important when dealing with overplotting.}

\item{exprs_use}{character, the name of the \code{assay} in the object whose values will
be plotted. Default is "logcounts".}

\item{knn_smooth}{logical, should the values be smoothed using a k-nearest neighbors
approach? Default is FALSE.}

\item{smoothing_k}{numeric, the number of neighbors to use when smoothing. Default is 10.}

\item{rng_seed}{numeric, the random number generator seed used when \code{plot_order = "random"}.}

\item{common_scale}{logical, should the points have a single, common color scale (TRUE) or 
should each panel have its own scale? Deafult is FALSE.}

\item{color_palette}{character, vector of colors to be interpolated across for the color
aesthetic. Default is NULL meaning a standard quantitative palette from \code{cellula}
will be used.}
}
\value{
a ggplot object showing a multi-panel plot of different features - one feature per panel -
   where there is either a common scale or an individual scale per panel.
}
\description{
Plots a multi-panel dimensionality reduction where points are colored 
by different features (one per panel)
}
\details{
This function is inspired by Seurat's \code{\link[Seurat]{FeaturePlot}}. 
    This function has two possible outcomes, controlled by the \code{common_scale} argument: 
    1) \code{common_scale = TRUE} returns a faceted plot where there is a single scale for color.
    This can be useful to compare gene expression across genes using color. 
    2) \code{common_scale = FALSE} returns a \code{patchwork} grid of feature plots, each with its own
    individual scale.
}
