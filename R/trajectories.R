#' Find trajectories
#'
#' Finds pseudo-temporal trajectories in a reduced dimensional space
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param dr character, the name of the \code{reducedDim} slot to use for trajectory
#'     estimation. Default is "PCA".
#' @param clusters character, the name of the \code{colData} column with cluster label
#'     indication.
#' @param method character, the method for estimation. Can be either \code{"slingshot"}
#'     or \code{"monocle"}. This will result in differences in the resulting object,
#'     see Details for more information.
#' @param ndims numeric, the number of dimensions to use in `dr`. If the number
#'     of columns in \code{dr} is < \code{ndims}, it will be used instead.
#' @param dr_embed character, the name of the \code{reducedDim} slot where curves should
#'     be embedded for plotting. If \code{method = "monocle"} and \code{dr_embed = "FR"}, 
#'     an alternative 2D embedding for the \code{monocle} trajectories is created (see
#'     Details).
#'     Default is NULL, meaning no embedding will be performed. 
#' @param start character, the name of the cluster to be used as starting point.
#'     Default is \code{"auto"}, implying an entropy-based method will be used to guess
#'     the best starting point.
#' @param Monocle_lg_control list or NULL (default). A list of control parameters
#'     for the \code{learn_graph()} function from \code{monocle3}. 
#'     See \code{?monocle3::learn_graph()}
#'     for more information. Only used when \code{method = "monocle"}.     
#' @param omega logical, should the \code{omega} method for MST calculation be used?
#'     Default is TRUE. See \code{?slingshot::getLineages} for more information.
#' @param omega_scale numeric, the value of the \code{omega_scale} parameter. 
#'     Default is 1.5. See \code{?slingshot::getLineages} for more information.          
#' @param do_de logical. Should differential expression across trajectories be
#'     performed? Default is FALSE.
#' @param batch_de character, the name of the \code{colData} column to be used as a
#'     blocking factor in the differential expression analysis. Default is NULL.
#' @param verbose logical, should progress messages be printed? Default is FALSE.    
#' @param BPPARAM a \code{BiocParallelParam} object. Default is NULL.
#'

#' @importFrom TSCAN perCellEntropy
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#'
#' @returns a SingleCellExperiment object with pseudotime results
#' 
#' @details
#' This function wraps calls to two popular trajectory estimation methods, i.e.
#' \code{slingshot} and \code{monocle3}. Most parameters are shared, and 
#' implementation is simplified for now (few parameters can be tweaked). In the 
#' future there will be more freedom to tweak parameters, although it bears
#' repeating that these are wrappers aimed at simplifying procedures.
#' 
#' For all methods the user is asked to define the dimensionality reduction slot 
#' via `dr`, and the number of dimensions via \code{ndims}. The user is also asked to
#' provide a starting cluster; if the choice is left to `auto`, the function will
#' use the entropy-based method as implemented in \code{TSCAN} to select maximum-
#' entropy cell clusters as starting points. Finally, both methods have the ability
#' to embed principal curves into a 2D representation of choice, albeit with 
#' slightly different results. 
#' 
#' The \code{slingshot} implementation allows the user to choose whether or not 
#' to use the \code{omega} method to separate disjointed trajectories, and to decide 
#' the \code{omega_scale}. Optionally, the user can run lineage-dependent differential
#' expression via \code{TSCAN::testPseudotime()} setting \code{do_de = TRUE}. 
#' 
#' The final result is a `SingleCellExperiment` object with 
#' some additional fields:
#' \itemize{
#'  \item{"slingPseudotime_N"}{ \code{colData} columns where N is any number >= 1. 
#'      These contain the pseudotemporal ordering of cells in a lineage, with NA
#'      being assigned to cells that do not belong to the lineage.}
#'  \item{"Slingshot_embedded_curves"}{ \code{metadata} list element containing segment 
#'      coordinates used to plot trajectories in 2D.}
#'  \item{"Slingshot_lineages"}{  \code{metadata} list element containing lineages (as 
#'      orderings of labels) used for metromap plotting.}   
#'  \item{"Slingshot_MST"}{ \code{metadata} list element containing the MST. Only
#'      added if \code{add_metadata = TRUE}.}
#'  \item{"Slingshot_curves"}{ \code{metadata} list element containing list of principal
#'      curve coordinates. Only added if `add_metadata` is `TRUE`.} 
#'  \item{"Slingshot_weights"}{ \code{metadata} list element containing a cell x lineage
#'      matrix of lineage-associated weights. Only added if `add_metadata` is 
#'      `TRUE`.}
#'  \item{"Slingshot_params"}{ \code{metadata} list element containing parameters for
#'      the \code{slingshot} call. Only added if `add_metadata` is `TRUE`.}
#'  \item{"pseudotime_DE"}{ \code{metadata} list element containing a list of DE results
#'      per lineage. Each result is a `DataFrame` object with `logFC`, `p.value`
#'      and `FDR` values for each gene. Only added if `do_de` is set to `TRUE}         
#' }
#' 
#' The \code{monocle3} implementation is rather simplified, with an important 
#' difference: it allows users to specify any reduced dimension rather than
#' just UMAP. This is in keeping with the evidence that UMAP reductions greatly
#' distorted. 
#'
#' When using PCA as a space for trajectory inference, the 2D 
#' embedding of the Monocle trajectories is re-calculated by picking the nearest
#' neighbors in PCA to the "waypoints" calculated by the algorithm. This can 
#' result in slightly more convoluted trajectories when visualized in UMAP, which
#' is attributable to the distortion.
#' 
#' To overcome this, an alternative embedding 
#' is available through \code{dr_embed = "FR"}. In this embedding the principal graph
#' of the trajectory is laid out using the Fruchterman-Reingold layout (hence
#' the name) and cells are randomly placed around their closest vertex (as 
#' calculated in PCA space) according to a 2D Gaussian distribution in which the
#' standard deviation is proportional to the square root of the number of cells
#' close to the vertex. Then, cell identities are ordered according to their
#' pseudotime value. Finally, this 2D embedding is given as an input to UMAP to 
#' optimize non-overlapping distribution of points. This results in a more pleasing
#' embedding in which cells are distributed along trajectories that were still
#' derived in high-dimensional space. This implementation was inspired by the
#' PAGA initialization for UMAP by Wolf and colleagues (2019).
#' 
#' The final result is a \code{SingleCellExperiment} object with 
#' some additional fields:
#' \itemize{
#'  \item{"monoclePseudotime"}{ \code{colData} column with a single pseudotime value
#'     for every cell.}
#'  \item{"Monocle_embedded_curves"}{ \code{metadata} list element containing segment 
#'      coordinates used to plot trajectories in 2D.}  
#'  \item{"Monocle_principal_graph"}{ \code{metadata} list element containing the 
#'      principal graph coordinates. Only added if `add_metadata` is set to 
#'      `TRUE`.} 
#'  \item{"Monocle_principal_graph_aux"}{ \code{metadata} list element containing 
#'      all the other objects used by \code{monocle3} for trajectory inference.
#'      Only added if \code{add_metadata = TRUE}}                
#' }
#' 
#' @author Giuseppe D'Agostino, and Stefan Boeing for the initial implementation of Monocle 3
#' 
#' @export

findTrajectories <- function(sce, dr = "PCA", clusters, method = "slingshot",
                             ndims = 20, dr_embed = NULL, start = "auto", 
                             Monocle_lg_control = NULL, omega = TRUE,
                             omega_scale = 1.5, invert = FALSE, do_de = FALSE, batch_de = NULL,
                             add_metadata = TRUE, verbose = FALSE, 
                             BPPARAM = SerialParam()) {

  if(!(start %in% colData(sce)[,clusters])) stop("Could not find the start cluster!")
    
  if(start == "auto") {
    if(!any(colnames(colData(sce)) == "entropy"))
      if(verbose) cat("Calculating per-cell entropy\n")
      sce$entropy = perCellEntropy(sce, BPPARAM)
    ent_means = lapply(split(colData(sce)$entropy, colData(sce)[,clusters]), mean)
    start = unique(colData(sce)[,clusters])[which.max(ent_means)]
  }
  
  if(ndims > ncol(reducedDims(sce, dr))) ndims = ncol(reducedDims(sce, dr))
  
  if(method == "slingshot"){
    
    sce = .getSlingshotTrajectories(sce = sce, dr = dr, ndims = ndims,
                                    clusters = clusters, start = start, 
                                    dr_embed = dr_embed, omega = omega, 
                                    omega_scale = omega_scale, do_de = do_de,
                                    batch_de = batch_de, verbose = verbose)
    
  } else if(method == "monocle") {
    
    sce = .getMonocleTrajectories(sce = sce, dr = dr, ndims = ndims,
                                  clusters = clusters, start = start, 
                                  Monocle_lg_control = Monocle_lg_control,
                                  dr_embed = dr_embed, invert = invert,
                                  add_metadata = add_metadata, verbose = verbose)
  }

  return(sce)
}


#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom monocle3 new_cell_data_set learn_graph principal_graph order_cells
#' @importFrom igraph V as_data_frame
#' @importFrom BiocNeighbors queryKNN
#' 
#' @noRd

.getMonocleTrajectories <- function(sce, 
                                    dr = "PCA", 
                                    ndims = 20, 
                                    clusters, 
                                    start,
                                    Monocle_lg_control = NULL,
                                    dr_embed = "UMAP",
                                    invert = FALSE,
                                    add_metadata = TRUE,
                                    verbose = TRUE){
  
  rd = rowData(sce)
  if(!("Symbol" %in% colnames(rd))) rd$Symbol = rownames(rd)
  rd$gene_short_name = rd$Symbol
  
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Creating CDS object"))
  
  cds <- new_cell_data_set(
    assay(sce, "counts"),
    cell_metadata = colData(sce),
    gene_metadata = rd
  )
  
  # Shoehorn any type of space into the UMAP slot
  reducedDim(cds, "UMAP") <- reducedDim(sce, dr)[,seq_len(ndims)]
  
  # only one partition - optional?
  recreate.partition <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partition) <- cds@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
  
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Adding cluster labels"))
  
  # Add cluster labels
  list_cluster = colData(sce)[, clusters, drop = TRUE]
  names(list_cluster) <- colnames(sce)
  cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  
  # Could be a space-holder, but essentially fills out louvain parameters
  
  cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Learning graph"))
  
  cds <- learn_graph(cds, use_partition = FALSE, 
                     learn_graph_control = Monocle_lg_control)
  
  # Set root cluster
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Finding start node"))
  
  root_cells <- as.vector(
    rownames(colData(sce))[which(colData(sce)[,clusters] == start)]
  )
  
  # Find starting vertex
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  rn = as.numeric(names(which.max(table(closest_vertex[root_cells,]))))
  root_pr_nodes <- V(principal_graph(cds)[["UMAP"]])$name[rn]
  
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Ordering cells"))
  
  cds <- order_cells( 
    cds, root_pr_nodes=root_pr_nodes)
  
  traj.coord <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  
  if (invert){
    PTmax <- max(traj.coord)
    traj.coord <- -1 * (traj.coord - PTmax)
    cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- traj.coord
  } 
  
  colData(sce)$monoclePseudotime = traj.coord
  
  # Embedding waypoint curves in 2D. We do this because we allow users to run 
  # Monocle3 in an arbitrary space that is not UMAP
  
  if(dr_embed == "FR") {
    if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Embedding in 2D"))
    
    frembed = .embedFR(cds = cds, sce = sce, dr = dr, ndims = ndims)
    
    metadata(sce)[["Monocle_embedded_curves"]] = frembed$segs
    
    #reducedDim(sce, "FR") = frembed$init
    
    reducedDim(sce, "UMAP_FR") = frembed$dr_embed
    
  } else if(!is.null(dr_embed) & (dr_embed != "FR")) {
    if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Embedding in 2D"))
    
    wp_coords = t(cds@principal_graph_aux$UMAP$dp_mst)
    sp_coords = reducedDim(sce, dr)[,seq_len(ndims)]
    
    nns = queryKNN(query = wp_coords,
                   X =  sp_coords,
                   k = 1)$index[,1]
    
    matched = data.frame("wp" = rownames(wp_coords), 
                         "sp" = rownames(sp_coords)[nns],
                         row.names = rownames(wp_coords))
    
    names(nns) = matched$sp
    wpsp_coords = reducedDim(sce, dr_embed)[nns,1:2]
    rownames(wpsp_coords) = rownames(wp_coords)
    
    segs = as_data_frame(cds@principal_graph$UMAP)
    
    segs$x0 = wpsp_coords[segs$from,1]
    segs$y0 = wpsp_coords[segs$from,2]
    segs$x1 = wpsp_coords[segs$to,1]
    segs$y1 = wpsp_coords[segs$to,2]
    
    segs$from = matched[segs$from, "sp"]
    segs$to = matched[segs$to, "sp"]
    
    metadata(sce)[["Monocle_embedded_curves"]] = segs
  }
  
  if(add_metadata) {
    if(verbose) message(paste0(blue("[TRAJ/Monocle3] "), "Adding metadata"))
    metadata(sce)[["Monocle_principal_graph"]] = cds@principal_graph$UMAP
    metadata(sce)[["Monocle_principal_graph_aux"]] = cds@principal_graph_aux$UMAP
  }
  
  if(verbose) message(paste0(blue("[TRAJ/Monocle3] "),"Done."))
  
  sce
}


#' @importFrom SummarizedExperiment colData 
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom BiocNeighbors queryKNN
#' @importFrom igraph V as_data_frame layout_with_fr subgraph graph_from_data_frame
#' @importFrom stats rnorm
#' 
#' @noRd

.embedFR <- function(cds, sce, dr, ndims) {
  
  wp_coords = t(cds@principal_graph_aux$UMAP$dp_mst)
  sp_coords = reducedDim(sce, dr)[,seq_len(ndims)]
  
  gv = cds@principal_graph$UMAP
  gv = as_data_frame(gv)
  cvert = as.data.frame(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex)
  verts = lapply(split(cvert, cvert$V1), rownames)
  names(verts) = names(V(cds@principal_graph$UMAP))
  gv$weight = (apply(gv[,1:2], 1, function(x) max(length(verts[[x[1]]]), length(verts[[x[2]]]))))
  l = layout_with_fr(graph_from_data_frame(gv))
  rownames(l) = names(V(graph_from_data_frame(gv)))
  
  coords = list()
  ptree =  cds@principal_graph_aux$UMAP$pr_graph_cell_proj_tree
  for(i in names(verts)) {
    xs = rnorm(n = length(verts[[i]]), mean = l[i,1], sd = 0.15*sqrt(length(verts[[i]])))
    ys = rnorm(n = length(verts[[i]]), mean = l[i,2], sd = 0.15*sqrt(length(verts[[i]])))
    pts = data.frame(x = xs, y = ys, row.names = verts[[i]]) #random assignment
    order_pt = names(V(subgraph(ptree, verts[[i]])))
    pts = pts[order_pt,]
    pts$order = order_pt
    coords[[i]] = pts
  }
  
  coords = do.call(rbind, coords)
  rownames(coords) = as.character(coords$order)
  init = as.matrix(coords[colnames(sce),1:2])
  
  frum = umap(init, 
              init = init, 
              n_neighbors = floor(sqrt(ncol(sce))), 
              min_dist = 1)
  
  segs = as_data_frame(cds@principal_graph$UMAP)
  
  nns = queryKNN(init, l, k = 1)
  matched = data.frame(wp = rownames(l), sp = nns$index, row.names = rownames(l))
  
  segs$x0 = frum[matched[segs$from,2],1]
  segs$y0 = frum[matched[segs$from,2],2]
  segs$x1 = frum[matched[segs$to,2],1]
  segs$y1 = frum[matched[segs$to,2],2]
  
  return(list(segs = segs, dr_embed = frum))
}

#' @importFrom slingshot slingshot slingLineages slingCurves embedCurves
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom TSCAN testPseudotime
#' 
#' @noRd

.getSlingshotTrajectories <- function(sce, dr = "PCA", 
                                      ndims = 20, 
                                      clusters, 
                                      dr_embed = NULL,
                                      start = "auto",
                                      omega = TRUE,
                                      omega_scale = 1.5,
                                      do_de = FALSE,
                                      batch_de = NULL,
                                      add_metadata = TRUE,
                                      verbose = FALSE,
                                      BPPARAM = SerialParam()){
  
  if(verbose) message(paste0(blue("[TRAJ/Slingshot] "),"Finding trajectories"))
  
  sce <- slingshot(sce,
                   reducedDim = dr,
                   clusterLabels = clusters,
                   start.clus = start,
                   omega = omega,
                   omega_scale = omega_scale)
  
  metadata(sce)[["Slingshot_lineages"]] = slingLineages(sce)
  
  if(!is.null(dr_embed)) {
    if(verbose) message(paste0(blue("[TRAJ/Slingshot] "),"Embedding curves"))
        ap = max(ncol(sce), 1000)
        embedded = lapply(slingCurves(embedCurves(sce, newDimRed = dr_embed, 
                                                  approx_points = ap)),
                          function(x) x$s)
        
        embedded = lapply(embedded, function(x) {
          data.frame(x0 = x[seq_len(ap-1),1], 
                     y0 = x[seq_len(ap-1),2], 
                     x1 = x[2:ap,1], 
                     y1 = x[2:ap,2])
                  })
        
        embedded = do.call(rbind, embedded)
        
        metadata(sce)[["Slingshot_embedded_curves"]] = embedded
  }
  
  if(do_de) {
    if(verbose) message(paste0(blue("[TRAJ/Slingshot] "),"Calculating DE along lineages"))
    if(!is.null(batch_de)) batch = factor(colData(sce)[,batch_de]) else batch = NULL
    
    sling_colnames = paste0("slingPseudotime_", seq_along(slingLineages(sce)))
    
    sling_tests <- lapply(sling_colnames, function(x)
      testPseudotime(assay(sce, "logcounts"),
                     pseudotime = colData(sce)[,x],
                     block = batch,
                     BPPARAM = BPPARAM)
    )
    
    names(sling_tests) <- names(slingLineages(sce))
    
    metadata(sce)[["pseudotime_DE"]] = sling_tests
  }
  
  if(add_metadata) {
    if(verbose) message(paste0(blue("[TRAJ/Slingshot] "),"Adding metadata"))
    metadata(sce)[["Slingshot_MST"]] = metadata(sce$slingshot)[["mst"]]
    metadata(sce)[["Slingshot_curves"]] = metadata(sce$slingshot)[["curves"]]
    metadata(sce)[["Slingshot_weights"]] = assay(sce$slingshot, "weights")
    metadata(sce)[["Slingshot_params"]] = metadata(sce$slingshot)[["slingParams"]]
  }
  
  sce$slingshot <- NULL
  
  sce
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
#' @param dr character, name of the `reducedDim` slot in which k-means clustering
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

makeMetacells  <- function(sce, w = 10, group = NULL, dr = "PCA", ndims = 20) {

  if(!is.null(group)) {
    gv = as.character(unique(colData(sce)[,group]))
    clustl = lapply(gv, function(x) {
      s = sce[,colData(sce)[,group] == x]
      spc = reducedDim(s, dr)[,seq_len(ndims)]
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
    dr = reducedDim(sce, dr)
    clust = kmeans(dr, centers = floor(ncol(sce)/w), iter.max = 50)
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


