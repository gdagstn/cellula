#' Find trajectories
#'
#' Finds pseudo-temporal trajectories in a reduced dimensional space
#'
#' @param sce a `SingleCellExperiment` object
#' @param space character, the name of the `reducedDim` slot to use for trajectory
#'     estimation. Default is "PCA".
#' @param clusters character, the name of the `colData` column with cluster label
#'     indication.
#' @param method character, the method for estimation. Currently only "slingshot".
#' @param ndims numeric, the number of dimensions to use in `space`. If the number
#'     of columns in `space` is < `ndims`, it will be used instead.
#' @param dr_embed character, the name of the `reducedDim` slot where curves should
#'     be embedded for plotting. Default is NULL, meaning no embedding will be
#'     performed.
#' @param start character, the name of the cluster to be used as starting point.
#'     Default is "auto", implying an entropy-based method will be used to guess
#'     the best starting point.
#' @param omega logical, should the `omega` method for MST calculation be used?
#'     Default is TRUE. See `?slingshot::getLineages` for more information.
#' @param omega_scale numeric, the value of the `omega_scale` parameter. 
#'     Default is 1.5. See `?slingshot::getLineages` for more information.          
#' @param do_de logical. Should differential expression across trajectories be
#'     performed? Default is FALSE.
#' @param batch_de character, the name of the `colData` column to be used as a
#'     blocking factor in the differential expression analysis. Default is NULL.
#' @param verbose logical, should progress messages be printed? Default is FALSE.    
#' @param BPPARAM a `BiocParallelParam` object. Default is NULL.
#'
#' @importFrom slingshot slingshot slingLineages slingCurves embedCurves
#' @importFrom TSCAN testPseudotime perCellEntropy
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#'
#' @returns a SingleCellExperiment object with slingshot outputs
#'
#' @export

findTrajectories <- function(sce, space = "PCA", clusters, method = "slingshot",
                             ndims = 20, dr_embed = NULL, start = "auto", omega = TRUE,
                             omega_scale = 1.5, do_de = FALSE, batch_de = NULL,
                             verbose = FALSE, BPPARAM = SerialParam()) {

  space = reducedDim(sce, space)[,seq_len(min(c(ndims, ncol(reducedDim(sce, space)))))]

  if(start == "auto") {
    if(!any(colnames(colData(sce)) == "entropy"))
      if(verbose) cat("Calculating per-cell entropy\n")
      sce$entropy = perCellEntropy(sce, BPPARAM)
    ent_means = lapply(split(colData(sce)$entropy, colData(sce)[,clusters]), mean)
    start = unique(colData(sce)[,clusters])[which.max(ent_means)]
  }

  if(method == "slingshot"){
    if(verbose) cat("Finding trajectories\n")
     sce <- slingshot(sce,
                      reducedDim = space,
                      clusterLabels = clusters,
                      start.clus = start,
                      omega = omega,
                      omega_scale = omega_scale)

     metadata(sce)[["Slingshot_lineages"]] = slingLineages(sce)

    if(!is.null(dr_embed)) {
      if(verbose) cat("Embedding curves\n")
      metadata(sce)[["Slingshot_embedded_curves"]] = lapply(slingCurves(embedCurves(sce, newDimRed = dr_embed)),
                                                            function(x) x$s)
     }

  if(do_de) {
    if(verbose) cat("Calculating DE along lineages\n")
    if(!is.null(batch_de)) batch = factor(colData(sce)[,batch_de]) else batch = NULL

    sling_colnames = paste0("slingPseudotime_", seq_along(slingLineages(sce)))
    
    sling_tests <- lapply(sling_colnames, function(x)
      testPseudotime(assay(sce, "logcounts"),
                            pseudotime = colData(sce)[,x],
                            block = batch,
                            BPPARAM = BPPARAM)
      )

      names(sling_tests) <- names(slingLineages(sce))

      metadata(sce)$pseudotime_DE = sling_tests
    }
  }

  return(sce)
}

#' Make metacells
#' 
#' Aggregate cells via k-means partitioning and summing
#' 
#' @param sce a SingleCellExperiment object
#' @param w numeric, the average number of cells that composes a metacell. Default
#'     is 10. 
#' @param group character, name of the `colData` column in which cells should be 
#'     grouped separately. Default is NULL.
#' @param space character, name of the `reducedDim` slot in which k-means clustering
#'     will be performed. Default is "PCA". 
#' @param ndims numeric, the number of dimensions of the space in which clustering
#'     will be performed. Default is 20.     
#'
#' @returns a SingleCellExperiment object with metacells and sample-level colData   
#' 
#' @importFrom stats kmeans
#' @importFrom SingleCellExperiment reducedDim SingleCellExperiment
#' @importFrom scuttle summarizeAssayByGroup
#' 
#' @export

makeMetacells  <- function(sce, w = 10, group = NULL, space = "PCA", ndims = 20) {

  if(!is.null(group)) {
    gv = as.character(unique(colData(sce)[,group]))
    clustl = lapply(gv, function(x) {
      s = sce[,colData(sce)[,group] == x]
      spc = reducedDim(s, space)[,seq_len(ndims)]
      clust = kmeans(spc, centers = floor(ncol(s)/w), iter.max = 50)
      memb = paste0(x, "_", clust$cluster)
      names(memb) = colnames(s)
      return(memb)
    })

    names(clustl) = gv
    memberships = unlist(clustl, use.names = TRUE)
    groups = rep(gv, lengths(clustl))
    names(memberships) = sapply(names(memberships), function(x) unlist(strsplit(x, split = "[.]"))[-1])
    ref = data.frame(groups = groups, membership = memberships, row.names = seq_along(memberships))
    ref = ref[!duplicated(ref),]
    rownames(ref) = ref$membership
  } else {
    space = reducedDim(sce, space)
    clust = kmeans(space, centers = floor(ncol(sce)/w), iter.max = 50)
    memberships = clust$cluster
    names(memberships) = colnames(sce)
  }

  agg = summarizeAssayByGroup(sce, ids = memberships[colnames(sce)], statistics = "sum")

  if(!is.null(group)) {
    colData(agg)[,group] = ref[colData(agg)[,"ids"], "groups"]
  }
  ret = SingleCellExperiment(assays = list(counts = assay(agg, "sum")),
                             colData = colData(agg),
                             rowData = rowData(agg))

  return(ret)
}


