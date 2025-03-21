% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doNormAndReduce.R
\name{paramSweepUMAP}
\alias{paramSweepUMAP}
\title{UMAP Parameter Sweep}
\usage{
paramSweepUMAP(
  sce,
  dr,
  ndims = 20,
  md = seq(0.1, 1, length.out = 10),
  nn = seq(5, 100, length.out = 20),
  seed = 11,
  consistent_seed = FALSE,
  verbose = FALSE,
  parallel_param = SerialParam()
)
}
\arguments{
\item{sce}{a \code{SingleCellExperiment} object}

\item{dr}{character, name of the \code{reducedDim} slot in \code{sce}}

\item{ndims}{numeric, number of dimensions (PCs) to be used for UMAP construction.
Default is 20.}

\item{md}{numeric vector, values for the minimum distance parameter grid}

\item{nn}{numeric vector, values for the number of neighbors parameter grid}

\item{seed}{numeric, the random number generator seed}

\item{consistent_seed}{logical, should the same seed be applied to all runs? 
Default is FALSE}

\item{verbose}{logical, display messages on progress? Default is \code{FALSE}.}

\item{parallel_param}{a \code{BiocParallel} object specifying the parallelization backend
to be used in some steps of the pipeline. Default is \code{SerialParam()},
meaning no parallelization will be used.}
}
\value{
a \code{list} containing two elements: 
    a \code{list} of UMAP embedding coordinates and a 
    \code{data.frame} with the parameter grid and its associated seed
}
\description{
Runs
}
