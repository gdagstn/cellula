#' SCE multi-resolution clustering
#'
#' Very simple wrapper to SNN graph and Louvain/Leiden clustering using
#' multiple resolutions
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param neighbors numeric, number of neighbors for SNN graph edge construction.
#'     Default is 10.
#' @param weighting_scheme character, the weighting scheme for SNN graph construction.
#'     One of \code{"jaccard"}, \code{"rank"}, \code{"number"}. Default is \code{"jaccard"}.
#' @param sweep_on character, the parameter used for sweeping. Can be \code{"clustering"},
#'     meaning values of \code{k} will be looped through as resolution, or \code{"SNN"},
#'     meaning values of \code{k} will be looped through as number of neighbors.
#' @param method character, the type of graph-based clustering to use. One of
#'     \code{"louvain"} or \code{"leiden"}. Default is \code{"louvain"}.
#' @param k numeric, vector of parameter sweep for graph construction or clustering.
#' @param space a matrix of lower dimensional embedding such as the PCA coordinates.
#'     if \code{NULL} (default), the \code{"PCA"} slot from \code{reducedDims(sce)}.
#' @param ndims numeric, the number of dimensions (columns of \code{space}) to use to
#'     build the SNN graph. Default is 20
#' @param calculate_modularity logical, should pairwise modularity between
#'     clusters be calculated? Default is \code{TRUE}
#' @param calculate_silhouette logical, should approximate silhouette widths be
#'     calculated? Default is \code{TRUE}
#' @param leiden_iterations numeric, the number of iterations of Leiden clustering.
#'     Default is 5.
#' @param save_graphs logical, should the \code{igraph} graph objects be saved in
#'     the metadata of \code{sce}? Default is FALSE.
#' @param prefix character, the prefix of the column names on \code{colData(sce)}
#'     where clustering results are stored. Default is \code{SNN_}.
#' @param verbose logical, should messages be written? Default is \code{FALSE}
#'
#' @return a \code{SingleCellExperiment} object with cluster memberships in the
#'     \code{colData} table, named according to prefix and respective value of \code{k}.
#'     Optionally, silhouette and/or modularity values are stored in the \code{metadata}
#'     slot of the \code{SingleCellExperiment} object, one for every value of \code{k}.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom crayon blue
#' @importFrom bluster makeSNNGraph pairwiseModularity approxSilhouette
#' @importFrom igraph cluster_louvain cluster_leiden
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#'
#' @export

makeGraphsAndClusters <- function(sce,
                                  neighbors = 10L,
                                  weighting_scheme = "jaccard",
                                  sweep_on = "clustering",
                                  method = "louvain",
                                  k = seq(0.1, 1, length.out = 6),
                                  space = NULL,
                                  ndims = 20L,
                                  calculate_modularity = TRUE,
                                  calculate_silhouette = TRUE,
                                  leiden_iterations = 5L,
                                  save_graphs = FALSE,
                                  prefix = "SNN_",
                                  verbose = FALSE) {
  
  #Sanity checks
  #Error prefix
  ep = "{cellula::makeGraphsAndClusters} - "
  
  if(!is(sce, "SingleCellExperiment")) 
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(!method %in% c("leiden", "louvain"))
    stop(paste0(ep, "method not recognized - must be one of \"leiden\" or \"louvain\""))
  if(!weighting_scheme %in%  c("jaccard", "rank", "number"))
    stop(paste0(ep, "method not recognized - must be one of \"jaccard\", \"rank\", or \"number\""))
  if(is.null(space)) {
    if(length(reducedDim(sce)) == 0) stop(paste0(ep, " there are no dimensionality reductions in the SingleCellExperiment object."))
    space = reducedDim(sce, "PCA")[,seq_len(ndims)]
  }

  # Case 1: parameter sweep on clustering resolution

  if(sweep_on == "clustering" & !is.null(neighbors)) {

    if(verbose) cat(blue("[CLU]"), "Creating SNN graph.\n")

    g = makeSNNGraph(space,
                     k = neighbors,
                     type = weighting_scheme)
    for(i in k){

      if(verbose) cat(blue("[CLU]"), "Clustering at resolution ", i, ".\n")

      if(method == "louvain") {
        cl = factor(cluster_louvain(g, resolution = i)$membership)
      } else if(method == "leiden") {
        cl = factor(cluster_leiden(g, objective_function = "modularity",
                                   n_iterations = leiden_iterations,
                                   resolution_parameter = i)$membership)
      }

      gname = paste0(prefix, i)
      colData(sce)[,gname] = cl

      if(verbose) cat(blue("[CLU]"), "Found", length(unique(cl)), "clusters.\n")

      if(calculate_modularity) {
        if(verbose) cat(blue("[CLU]"), "Calculating pairwise modularity.\n")
        metadata(sce)[[paste0("modularity_", gname)]] = pairwiseModularity(g, clusters = cl, as.ratio = TRUE)
      }

      if(calculate_silhouette) {
        if(verbose) cat(blue("[CLU]"), "Calculating approximate silhouette widths.\n")
        silhouette <- as.data.frame(approxSilhouette(space, clusters = cl))
        silhouette$closest <- factor(ifelse(silhouette$width > 0, cl, silhouette$other))
        silhouette$cluster <- cl
        metadata(sce)[[paste0("silhouette_", gname)]] = silhouette
      }
      
      if(save_graphs) {
        metadata(sce)[[paste0("SNN_", neighbors, "_", weighting_scheme)]] = g
      }
    }
    # Case 2: parameter sweep on SNN neighbor number
  } else if(sweep_on == "SNN" & !is.null(k)) {

    k = floor(k)

    for(i in k) {
      if(verbose) cat(blue("[CLU]"), "Creating SNN graph with k =", i, "neighbors.\n")
      g = makeSNNGraph(space,
                       k = i,
                       type = weighting_scheme)
      
      if(save_graphs) {
        metadata(sce)[[paste0("SNN_", i, "_", weighting_scheme)]] = g
      }
      
      if(method == "louvain") {
        cl = factor(cluster_louvain(g)$membership)
      } else if(method == "leiden") {
        cl = factor(cluster_leiden(g, objective_function = "modularity",
                                   n_iterations = leiden_iterations)$membership)
      }
      gname = paste0(prefix, i)
      colData(sce)[,gname] = cl
      if(verbose) cat("Found", length(unique(cl)), "clusters.\n")

      if(calculate_modularity) {
        if(verbose) cat(blue("[CLU]"), "Calculating pairwise modularity.\n")
        metadata(sce)[[paste0("modularity_", gname)]] = pairwiseModularity(g, clusters = cl, as.ratio = TRUE)
      }

      if(calculate_silhouette) {
        if(verbose) cat(blue("[CLU]"), "Calculating approximate silhouette widths.\n")
        sil <- as.data.frame(approxSilhouette(space, clusters = cl))
        sil$closest <- factor(ifelse(sil$width > 0, cl, sil$other))
        sil$cluster <- cl
        metadata(sce)[[paste0("silhouette_", gname)]] = sil
      }
    }
  }

  return(sce)
}


#' SCE meta-clustering
#'
#' Applies a meta-clustering procedure through cluster linking to identify
#' consensus clusters
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param clusters character, the column names where cluster assignments can be
#'    found in \code{colData(sce)}
#' @param threshold numeric, the score threshold to determine metacluster assignment.
#'    Default is 0.5, must be between 0 and 1.
#' @param do_plot logical, should the metacluster plot be printed? Default is \code{TRUE}.
#' @param denominator character, one of \code{"min"}, \code{"max"}, \code{"union"}. 
#'    Default is \code{"union".} See \code{?bluster::linkClusters} for details.
#'
#' @return a SingleCellExperiment object with three additional columns:
#'    - metacluster_max indicating the most frequent metacluster assignment
#'    - metacluster_score for the assignment score (between 0 and 1)
#'    - metacluster_ok to test whether the score is above threshold
#'
#' @importFrom bluster linkClusters
#' @importFrom igraph cluster_louvain E groups
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export

metaCluster <- function(sce,
                        clusters,
                        threshold = 0.5,
                        do_plot = TRUE,
                        denominator = "union") {

  #Sanity checks
  #Error prefix
  ep = "{cellula::metaCluster} - "
  
  if(threshold > 1 | threshold < 0) stop(paste0(ep, "threshold must be between 0 and 1"))
  if(any(!clusters %in% colnames(colData(sce)))) stop(paste0(ep, "some cluster column names were not found in colData"))
  if(!is(sce, "SingleCellExperiment")) stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if((!denominator %in% c("union", "max", "min"))) stop(paste0(ep, "denominator must be one of \"union\", \"max\", \"min\""))

  linked <- linkClusters(colData(sce)[,clusters], denominator = denominator)
  meta <- cluster_louvain(linked)

  if(do_plot) {
    pal = rainbow(length(groups(meta)), alpha = 0.3)
    plot(linked, vertex.color = pal[as.numeric(as.factor(meta$membership))],
       edge.width = E(linked)$weight * 2, vertex.size = 5, vertex.label.cex = 0.4)
    legend("topleft", pt.bg = pal, pch = 21, legend = 1:10, bty = "n", title = "Metaclusters")
  }

  meta_memberships = data.frame(row.names = meta$names, "metacluster" = meta$membership)

  clusterlabels = lapply(clusters, function(x) paste0(x, ".", colData(sce)[,x]))

  metalabels = lapply(clusterlabels, function(x) meta_memberships[x,])
  metalabeldf = as.data.frame(do.call(cbind, metalabels))
  colnames(metalabeldf) = clusters
  metafuzzy = apply(metalabeldf, 1, function(x) max(table(as.numeric(x)))/ncol(metalabeldf))

  sce$metacluster_max = factor(apply(metalabeldf, 1, function(x) names(table(x))[which.max(table(x))]))
  sce$metacluster_score = metafuzzy
  sce$metacluster_ok = metafuzzy < threshold
  metadata(sce)$metaclusters <- metalabeldf

  return(sce)
}

