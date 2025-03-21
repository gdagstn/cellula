% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doQC.R
\name{doQC}
\alias{doQC}
\title{SCE QC sub-pipeline}
\usage{
doQC(
  sce,
  batch = NULL,
  name,
  discard = TRUE,
  nmad = 3,
  subset_mito = TRUE,
  subset_ribo = TRUE,
  subset_malat1 = TRUE,
  detect_doublets = TRUE,
  run_emptydrops = FALSE,
  emptydrops_cutoff = "auto",
  emptydrops_alpha = 0.01,
  verbose = TRUE,
  save_plots = TRUE,
  path = "./",
  parallel_param = SerialParam()
)
}
\arguments{
\item{sce}{a \code{SingleCellExperiment} object}

\item{batch}{character, the name of the column in \code{colData(sce)} with batch labels.
Default is \code{NULL} meaning no batches will be considered, and data will be
processed as a single batch.}

\item{name}{character, the name of the file/folder.}

\item{discard}{logical, should values that do not meet QC thresholds be discarded?
Default is \code{TRUE}.}

\item{nmad}{numeric, the number of median absolute deviations to use for
outlier detection. Default is \code{3}.}

\item{subset_mito}{logical, should mitochondrial transcripts be used for QC?
Default is \code{TRUE}}

\item{subset_ribo}{logical, should ribosomal transcripts be used for QC?
Default is \code{TRUE}}

\item{subset_malat1}{logical, should MALAT1 transcripts be used for QC?
Default is \code{TRUE}}

\item{detect_doublets}{logical, should \code{scDblFinder} be run? Default is \code{TRUE}}

\item{run_emptydrops}{logical, should `emptyDrops` be run? Default is \code{TRUE}}

\item{emptydrops_cutoff}{either \code{"auto"} (default, barcode rank inflection point)
or a numeric. Cells with total reads below this cutoff are used to calculate
ambient RNA profiles and are removed.}

\item{emptydrops_alpha}{numeric, the FDR threshold to call an empty barcode.}

\item{verbose}{logical, display messages on progress? Default is \code{FALSE}}

\item{save_plots}{logical, should plots be drawn and saved? Default is \code{TRUE}}

\item{path}{character, the path to save the plots. Default is the local working directory ("./")}

\item{parallel_param}{a \code{BiocParallel} object specifying the parallelization backend
to be used in some steps of the pipeline. Default is \code{SerialParam()}
meaning no parallelization will be used.}
}
\value{
a \code{SingleCellExperiment} object with QC including doublet assignment
}
\description{
Pipeline for automatic processing and integration of SingleCellExperiment objects
}
