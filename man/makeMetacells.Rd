% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectories.R
\name{makeMetacells}
\alias{makeMetacells}
\title{Make metacells}
\usage{
makeMetacells(sce, w = 10, group = NULL, dr = "PCA", ndims = 20)
}
\arguments{
\item{sce}{a SingleCellExperiment object}

\item{w}{numeric, the average number of cells that composes a metacell. Default
is 10.}

\item{group}{character, name of the `colData` column in which cells should be 
grouped separately. Default is NULL.}

\item{dr}{character, name of the `reducedDim` slot in which k-means clustering
will be performed. Default is "PCA".}

\item{ndims}{numeric, the number of dimensions of the space in which clustering
will be performed. Default is 20.}
}
\value{
a SingleCellExperiment object with metacells and sample-level colData
}
\description{
Aggregate cells via k-means partitioning and summing
}
