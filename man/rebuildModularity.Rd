% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeGraphsAndClusters.R
\name{rebuildModularity}
\alias{rebuildModularity}
\title{Rebuild modularity}
\usage{
rebuildModularity(
  sce,
  dr = "PCA",
  neighbors = NULL,
  ndims = 20,
  labels,
  weighting_scheme = "jaccard"
)
}
\arguments{
\item{sce}{a \code{SingleCellExperiment} object}

\item{dr}{character, the name of the \code{reducedDim} slot to be used.
Default is \code{"PCA"}}

\item{neighbors}{numeric, the number of neighbors to be used for SNN graph 
construction. Default is \code{NULL} meaning the square root of the 
number of cells will be used.}

\item{ndims, }{numeric, the number of dimensions to use from \code{dr}. Default is 20.}

\item{labels}{character, the name of the \code{colData} column with cluster
labels to be used.}

\item{weighting_scheme}{character, one of \code{"jaccard"} (default), 
\code{rank} or \code{"min"}}
}
\value{
a \code{SingleCellExperiment} object with a modularity matrix in the
    \code{metadata} slot, named \code{modularity_$labels} where $labels is the
    user-defined column name.
}
\description{
Rebuilds a pairwise modularity matrix using user-defined labels
}
\details{
This function can be used to approximate a pairwise modularity matrix
    in case labels are not generated through the \code{\link{makeGraphsAndClusters}}
    function. It is considered to be _approximate_ because it requires the 
    user to define some parameters (neighbors, dimensions, weighting scheme)
    to rebuild a Shared Nearest Neighbor graph. If the user has defined labels
    without clustering (e.g. by running \code{\link{assignIdentities}}) there is no
    guarantee that these parameters are the ones that result in a clustering
    that easily translate to the same labels. This is a convenience function.
    The output can be used by the \code{"MODDPT"} method in the trajectory
    inference module in case there are no previously available modularity 
    matrices.
}
