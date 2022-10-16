

#' SCE multi-resolution clustering
#'
#' Very simple wrapper to SNN graph and Louvain/Leiden clustering using
#' multiple resolutions
#'
#' @param sce a SingleCellExperiment object
#' @param neighbors numeric, number of neighbors for SNN graph edge construction.
#'     Default is 10.
#' @param weighting_scheme character, the weighting scheme for SNN graph construction.
#'     One of "jaccard", "rank", "number". Default is "jaccard".
#' @param sweep_on character, the parameter used for sweeping. Can be "clustering",
#'     meaning values of `k` will be looped through as resolution, or "SNN",
#'     meaning values of `k`will be looped through as number of neighbors.
#' @param method character, the type of graph-based clustering to use. One of
#'     "louvain" or leiden". Default is "louvain".
#' @param k numeric, vector of parameter sweep for graph cosntruction or clustering.
#' @param space a matrix of lower dimensional embedding such as the PCA coordinates.
#'     if NULL (default), the "PCA" slot from `reducedDims(sce)`
#' @param ndims numeric, the number of dimensions (columns of `space`) to use to
#'     build the SNN graph. Default is 20
#' @param calculate_modularity logical, should pairwise modularity between
#'     clusters be calculated? Default is TRUE
#' @param calculate_silhouette logical, should approximate silhouette widths be
#'     calculated? Default is TRUE
#' @param leiden_iterations numeric, the number of iterations of Leiden clustering.
#'     Default is 5.
#' @param prefix character, the prefix of the column names on `colData`
#'     where clustering results are stored. Default is "SNN_".
#' @param verbose logical, should messages be written? Default is FALSE
#' @param BPPARAM a `BiocParallel` parameter object for graph construction
#'
#' @return a `SingleCellExperiment` object with cluster memberships in the
#'     `colData` table, named according to prefix and respective value of `k`.
#'     Optionally, silhouette and/or modularity values are stored in the `metadata`
#'     slot of the `SingleCellExperiment` object, one for every value of `k`.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom crayon blue
#' @importFrom bluster makeSNNGraph pairwiseModularity approxSilhouette
#' @importFrom igraph cluster_louvain cluster_leiden
#' @importFrom SingleCellExperiment reducedDim reducedDim<- counts logcounts
#' @importFrom BiocParallel SerialParam
#'
#'
#' @export

makeGraphsAndClusters <- function(sce,
                                  neighbors = 10,
                                  weighting_scheme = "jaccard",
                                  sweep_on = "clustering",
                                  method = "louvain",
                                  k = seq(0.1, 1, length.out = 6),
                                  space = NULL,
                                  ndims = 20,
                                  calculate_modularity = TRUE,
                                  calculate_silhouette = TRUE,
                                  leiden_iterations = 5,
                                  prefix = "SNN_",
                                  verbose = FALSE,
                                  BPPARAM = SerialParam()) {

  if(is.null(space)) space = reducedDim(sce, "PCA")[,seq_len(ndims)]


  # Case 1: parameter sweep on clustering resolution

  if(sweep_on == "clustering" & !is.null(neighbors)) {

    if(verbose) cat(blue("[CLU]"), "Creating SNN graph.\n")

    g = makeSNNGraph(space,
                     k = neighbors,
                     type = weighting_scheme,
                     BPPARAM = BPPARAM)
    for(i in k){

      if(verbose) cat(blue("[CLU]"), "Clustering at resolution ", i, ".\n")

      if(method == "louvain") {
        cl = factor(cluster_louvain(g, resolution = i)$membership)
      } else if(method == "leiden") {
        cl = factor(cluster_leiden(g, objective_function = "CPM",
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
    }
    # Case 2: parameter sweep on SNN neighbor number
  } else if(sweep_on == "SNN" & !is.null(k)) {

    k = floor(k)

    for(i in k) {
      if(verbose) cat(blue("[CLU]"), "Creating SNN graph with k =", i, "neighbors.\n")
      g = makeSNNGraph(space,
                       k = i,
                       type = weighting_scheme,
                       BPPARAM = BPPARAM)
      if(method == "louvain") {
        cl = factor(cluster_louvain(g)$membership)
      } else if(method == "leiden") {
        cl = factor(cluster_leiden(g, objective_function = "CPM",
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
