% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeGraphsAndClusters.R
\name{metaCluster}
\alias{metaCluster}
\title{SCE meta-clustering}
\usage{
metaCluster(
  sce,
  clusters,
  threshold = 0.5,
  do_plot = TRUE,
  denominator = "union"
)
}
\arguments{
\item{sce}{a \code{SingleCellExperiment} object}

\item{clusters}{character, the column names where cluster assignments can be
found in \code{colData(sce)}}

\item{threshold}{numeric, the score threshold to determine metacluster assignment.
Default is 0.5, must be between 0 and 1.}

\item{do_plot}{logical, should the metacluster plot be printed? Default is \code{TRUE}.}

\item{denominator}{character, one of \code{"min"}, \code{"max"}, \code{"union"}. 
Default is \code{"union".} See \code{\link[bluster]{linkClusters}} for details.}
}
\value{
a SingleCellExperiment object with three additional columns:
   - metacluster_max indicating the most frequent metacluster assignment
   - metacluster_score for the assignment score (between 0 and 1)
   - metacluster_ok to test whether the score is above threshold
}
\description{
Applies a meta-clustering procedure through cluster linking to identify
consensus clusters
}
