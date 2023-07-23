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
#'     See \code{?monocle3::learn_graph}
#'     for more information. Only used when \code{method = "monocle"}.     
#' @param omega logical, should the \code{omega} method for MST calculation be used?
#'     Default is TRUE. See \code{?slingshot::getLineages} for more information.
#' @param omega_scale numeric, the value of the \code{omega_scale} parameter. 
#'     Default is 1.5. See \code{?slingshot::getLineages} for more information.  
#' @param invert logical. Should the pseuodtime vector be inverted? Only valid
#'     for monocle3. Default is FALSE.            
#' @param do_de logical. Should differential expression across trajectories be
#'     performed? Default is FALSE.
#' @param batch_de character, the name of the \code{colData} column to be used as a
#'     blocking factor in the differential expression analysis. Default is NULL.
#' @param add_metadata logical, should additional data from trajectory inference be
#'     added to the \code{metadata(sce)}? Default is TRUE.
#' @param verbose logical, should progress messages be printed? Default is FALSE.    
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
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
#' use an entropy-based method to select maximum-entropy cell clusters as starting 
#' points. Finally, both methods have the ability to embed principal curves into 
#' a 2D representation of choice, albeit with different results. 
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
                             add_metadata = TRUE, verbose = FALSE) {
  
  ## Sanity checks
  # Error prefix
  ep = "{cellula::findTrajectories()} - "
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(ep, "Must provide a SingleCellExperiment object"))
  if(!(method %in% c("slingshot", "monocle", "PAGADPT"))) 
    stop(paste0(ep, "method not recognized - must be one of \"slingshot\", \"monocle\", or \"PAGADPT\""))
  if(!(clusters %in% colnames(colData(sce)))) 
    stop(paste0(ep,"clusters column not found in the colData of the object"))
  if(!is(colData(sce)[,clusters], "character") & !is(colData(sce)[,clusters], "factor")) 
    stop(paste0(ep,"clusters column should contain a factor or a character"))
  if(!(dr %in% reducedDimNames(sce)))
    stop(paste0(ep,"dr reduction not found among the reducedDims of the object"))
  if(!is.null(dr_embed)){ 
    if(dr_embed != "FR" & !(dr_embed %in% reducedDimNames(sce))) 
      stop(paste0(ep,"dr_embed reduction not found among the reducedDims of the object"))
  }
  if(start != "auto" & !(start %in% colData(sce)[,clusters])) 
    stop(paste0(ep,"Could not find the start cluster in clusters"))
  
  # Start parameter logging - not fully implemented
  # TO DO
  # --------------------------------------------- #
  
  if(start == "auto") {
    if(!any(colnames(colData(sce)) == "entropy"))
      if(verbose) cat("Calculating per-cell entropy\n")
      sce$entropy = apply(counts(sce), 2, .getEntropy) 
    ent_means = lapply(split(colData(sce)$entropy, colData(sce)[,clusters]), mean)
    start = unique(colData(sce)[,clusters])[which.max(ent_means)]
  }
  
  if(ndims > ncol(reducedDim(sce, dr))) ndims = ncol(reducedDim(sce, dr))
  
  # Trajectory inference module
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
    
  } else if(method == "PAGADPT") { # UNDOCUMENTED (WIP)
    sce = .getDPTtrajectories(sce = sce, dr = dr, ndims = ndims, 
                              clusters = clusters, start = start, 
                              add_metadata = add_metadata, verbose = verbose)
  }

  return(sce)
}


#' @importFrom SummarizedExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom igraph V as_data_frame
#' @importFrom BiocNeighbors queryKNN
#' @importFrom S4Vectors metadata metadata<-
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
  
  ep = "{cellula::.getMonocleTrajectories()} - "
  
  if(!"monocle3" %in% rownames(installed.packages()))
    stop(paste0(.redm(ep), "the `monocle3` package must be installed first.\n
                Run `BiocManager::install(\"cole-trapnell-lab/monocle3\") to use this function."))
  
  rd = rowData(sce)
  if(!("Symbol" %in% colnames(rd))) rd$Symbol = rownames(rd)
  rd$gene_short_name = rd$Symbol
  
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Creating CDS object"))
  
  cds <- monocle3::new_cell_data_set(
                    assay(sce, "counts"),
                    cell_metadata = colData(sce),
                    gene_metadata = rd
                  )
                  
  # Shoehorn any type of space into the UMAP slot
  reducedDim(cds, "UMAP") <- reducedDim(sce, dr)[,seq_len(ndims)]
  
  # only one partition - NEED TO REVISIT
  recreate.partition <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partition) <- cds@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
  
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Adding cluster labels"))
  
  # Add cluster labels
  list_cluster = colData(sce)[, clusters, drop = TRUE]
  names(list_cluster) <- colnames(sce)
  cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
  
  # Could be a space-holder, but essentially fills out louvain parameters
  cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Learning graph"))
  
  cds <- monocle3::learn_graph(cds, use_partition = FALSE, 
                               learn_graph_control = Monocle_lg_control)
  
  # Set root cluster
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Finding start node"))
  
  root_cells <- as.vector(
    rownames(colData(sce))[which(colData(sce)[,clusters] == start)]
  )
  
  # Find starting vertex
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  rn = as.numeric(names(which.max(table(closest_vertex[root_cells,]))))
  root_pr_nodes <- V(monocle3::principal_graph(cds)[["UMAP"]])$name[rn]
  
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Ordering cells"))
  
  cds <- monocle3::order_cells(cds, root_pr_nodes=root_pr_nodes)
  
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
    if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Embedding in 2D"))
    
    frembed = .embedFR(cds = cds, sce = sce, dr = dr, ndims = ndims)
    
    metadata(sce)[["Monocle_embedded_curves"]] = frembed$segs
    
    #reducedDim(sce, "FR") = frembed$init
    
    reducedDim(sce, "UMAP_FR") = frembed$dr_embed
    
  } else if(!is.null(dr_embed) & (dr_embed != "FR")) {
    if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Embedding in 2D"))
    
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
    if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "), "Adding metadata"))
    metadata(sce)[["Monocle_principal_graph"]] = cds@principal_graph$UMAP
    metadata(sce)[["Monocle_principal_graph_aux"]] = cds@principal_graph_aux$UMAP
  }
  
  if(verbose) message(paste0(.bluem("[TRAJ/Monocle3] "),"Done."))
  
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
  
  # Vertices of the principal graph in UMAP
  gv = as_data_frame(cds@principal_graph$UMAP)
  cvert = as.data.frame(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex)
  verts = lapply(split(cvert, cvert$V1), rownames)
  names(verts) = names(V(cds@principal_graph$UMAP))
  
  # Weight by the number of cells closest to each node - more cells, closer points
  gv$weight = apply(gv[,1:2], 1, function(x) max(length(verts[[x[1]]]), length(verts[[x[2]]])))
  
  # Fruchterman-Rheingold layout
  l = layout_with_fr(graph_from_data_frame(gv))
  rownames(l) = names(V(graph_from_data_frame(gv)))
  
  coords = list()
  ptree =  cds@principal_graph_aux$UMAP$pr_graph_cell_proj_tree
  
  traj.coord <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  
  # Randomize position of cells on the FR layout but keep order consistent with PT
  for(i in names(verts)) {
    xs = rnorm(n = length(verts[[i]]), mean = l[i,1], sd = 0.15*sqrt(length(verts[[i]])))
    ys = rnorm(n = length(verts[[i]]), mean = l[i,2], sd = 0.15*sqrt(length(verts[[i]])))
    pts = data.frame(x = xs, y = ys, row.names = verts[[i]]) #random assignment
    
    order_pt = names(V(subgraph(ptree, verts[[i]])))
    
    extreme_pts = traj.coord[c(order_pt[1], order_pt[length(order_pt)])]
    if(diff(extreme_pts) < 0) order_pt = rev(order_pt)
      
    pts = pts[order_pt,]
    pts$order = order_pt
    coords[[i]] = pts
  }
  
  coords = do.call(rbind, coords)
  rownames(coords) = as.character(coords$order)
  init = as.matrix(coords[colnames(sce),1:2])
  
  # Reassign points using UMAP
  frum = umap(init, 
              init = init, 
              n_neighbors = floor(sqrt(ncol(sce))), 
              min_dist = 1)
  
  segs = as_data_frame(cds@principal_graph$UMAP)
  
  # Segments of the principal graph on UMAP
  nns = queryKNN(init, l, k = 1)
  matched = data.frame(wp = rownames(l), sp = nns$index, row.names = rownames(l))
  
  segs$x0 = frum[matched[segs$from,2],1]
  segs$y0 = frum[matched[segs$from,2],2]
  segs$x1 = frum[matched[segs$to,2],1]
  segs$y1 = frum[matched[segs$to,2],2]
  
  return(list(segs = segs, dr_embed = frum))
}

#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata metadata<-
#'
#' @noRd

.getSlingshotTrajectories <- function(sce, dr = "PCA", 
                                      ndims = 20, 
                                      clusters, 
                                      dr_embed = NULL,
                                      start,
                                      omega = TRUE,
                                      omega_scale = 1.5,
                                      do_de = FALSE,
                                      batch_de = NULL,
                                      add_metadata = TRUE,
                                      verbose = FALSE,
                                      BPPARAM = SerialParam()){
  
  ep = "{cellula::.getSlingshotTrajectories()} - "
  
  if(!"slingshot" %in% rownames(installed.packages()))
    stop(paste0(ep, "the `slingshot` package must be installed first.\n
                Run `BiocManager::install(\"slingshot\") to use this function."))
  
  if(verbose) message(paste0(.bluem("[TRAJ/Slingshot] "),"Finding trajectories"))
  
  sce <- slingshot::slingshot(sce,
                             reducedDim = dr,
                             clusterLabels = clusters,
                             start.clus = start,
                             omega = omega,
                             omega_scale = omega_scale)
  
  metadata(sce)[["Slingshot_lineages"]] = slingshot::slingLineages(sce)
  
  if(!is.null(dr_embed)) {
    if(verbose) message(paste0(.bluem("[TRAJ/Slingshot] "),"Embedding curves"))
        ap = max(ncol(sce), 1000)
        embedded = lapply(slingshot::slingCurves(
          slingshot::embedCurves(sce, newDimRed = dr_embed, approx_points = ap)),
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
    if(verbose) message(paste0(.bluem("[TRAJ/Slingshot] "),"Calculating DE along lineages"))
    if(!is.null(batch_de)) batch = factor(colData(sce)[,batch_de]) else batch = NULL
    
    sling_colnames = paste0("slingPseudotime_", seq_along(slingshot::slingLineages(sce)))
    
    sling_tests <- lapply(sling_colnames, function(x)
      TSCAN::testPseudotime(assay(sce, "logcounts"),
                     pseudotime = colData(sce)[,x],
                     block = batch,
                     BPPARAM = BPPARAM)
    )
    
    names(sling_tests) <- names(slingshot::slingLineages(sce))
    
    metadata(sce)[["pseudotime_DE"]] = sling_tests
  }
  
  if(add_metadata) {
    if(verbose) message(paste0(.bluem("[TRAJ/Slingshot] "),"Adding metadata"))
    metadata(sce)[["Slingshot_MST"]] = metadata(sce$slingshot)[["mst"]]
    metadata(sce)[["Slingshot_curves"]] = metadata(sce$slingshot)[["curves"]]
    metadata(sce)[["Slingshot_weights"]] = assay(sce$slingshot, "weights")
    metadata(sce)[["Slingshot_params"]] = metadata(sce$slingshot)[["slingParams"]]
  }
  
  sce$slingshot <- NULL
  
  sce
}

#' @importFrom igraph graph_from_adjacency_matrix E V mst degree all_shortest_paths
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData 
#' @importFrom S4Vectors metadata metadata<-
#' 
#' @noRd

.getDPTtrajectories <- function(sce, 
                                dr = "PCA", 
                                ndims = 20, 
                                clusters, 
                                #dr_embed = NULL, 
                                start, 
                                add_metadata = TRUE,
                                verbose = FALSE){#,
                                #BPPARAM = SerialParam()) {
  
  ep = "{cellula::.getDPTtrajectories()} - "
  
  if(!"destiny" %in% rownames(installed.packages()))
    stop(paste0(.redm(ep), "the `destiny` package must be installed first.\n
                Run `BiocManager::install(\"destiny\") to use this function."))
  
  if(verbose) message(paste0(.bluem("[TRAJ/PAGA-DPT] "),"Building MST on modularity graph"))
  
  modslot = paste0("modularity_", clusters)
  
  if(!modslot %in% names(metadata(sce)))
    stop(paste0(ep, "Modularity matrix not found! Pairwise modularity should have been 
         calculated in the clustering step and saved to the metadata."))
  
  modmat = metadata(sce)[[modslot]]
  mgr <- graph_from_adjacency_matrix(modmat, mode = "upper", 
                                     weighted = TRUE, diag = FALSE)
  msgr = mst(mgr, weights = 1/E(mgr)$weight)
  outd = degree(msgr, v = V(msgr), mode = "out")
  outnodes = names(outd)[outd == 1]
  paths = all_shortest_paths(msgr, from = start, to = outnodes)
  paths = lapply(paths$res, function(x) as.numeric(x))
  
  dptlist = dmlist = list()
  
  if(verbose) message(paste0(.bluem("[TRAJ/PAGA-DPT] "),"Calculating diffusion maps and DPT"))
  
  for(i in seq_along(paths)) {
    if(verbose) {
      path_text = paste(paths[[i]], collapse = "-->")
    }
    curr = sce[,colData(sce)[,clusters] %in% paths[[i]]]
    reducedDim(curr, "pca") = reducedDim(curr, dr)
    dmlist[[i]] = destiny::DiffusionMap(data = reducedDim(curr, "pca"), verbose = verbose)
    dptlist[[i]] = destiny::DPT(dmlist[[i]])$dpt
    names(dptlist[[i]]) = colnames(curr)
  }
  
  for(i in seq_along(paths)) {
    colData(sce)[[paste0("DPTpseudotime_", i)]] = NA
    colData(sce)[names(dptlist[[i]]), paste0("DPTpseudotime_", i)] = dptlist[[i]]
  }
  
  if(add_metadata) {
    if(verbose) message(paste0(.bluem("[TRAJ/PAGA-DPT] "),"Adding metadata"))
    metadata(sce)[["PAGADPT_MST"]] = msgr
    metadata(sce)[["PAGADPT_lineages"]] = paths
    metadata(sce)[["PAGADPT_DiffusionMaps"]] = dmlist
  }
  
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
 
  ep = "{cellula::makeMetacells()} - "
  
  if(!is(sce, "SingleCellExperiment"))
    stop(paste0(.redm(ep), "Must provide a SingleCellExperiment object"))
    if(!(group %in% colnames(colData(sce)))) 
    stop(paste0(.redm(ep),"group column not found in the colData of the object"))
  if(!is(colData(sce)[,group], "character") & !is(colData(sce)[,group], "factor")) 
    stop(paste0(.redm(ep),"group column should contain a factor or a character"))
  if(!(dr %in% reducedDimNames(sce)))
    stop(paste0(.redm(ep),"dr reduction not found among the reducedDims of the object"))
  if(ndims > ncol(reducedDim(sce, dr)))
    stop(paste0(.redm(ep), "ndims is more than the available dimensions for the selected dr slot."))
  
  # Start parameter logging - not fully implemented
  # TO DO
  # --------------------------------------------- #
  
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


#' @noRd
.getEntropy <- function(x) {
    x <- x[x > 0]
    x <- x/sum(x)
    -sum(x * log(x))
}
