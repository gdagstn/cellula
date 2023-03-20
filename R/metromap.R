#' Plot Metro Map
#'
#' Plots a metro map-like representation of trajectories
#'
#' @param sce a `SingleCellExperiment` class object
#' @param dimred character, the name of the 2D dimensional reduction to use as
#'     base for plotting
#' @param snap logical, should points and lines be snapped to a regular grid?
#'     default is TRUE.
#' @param res numeric, the resolution of the grid for snapping. Higher = smaller
#'     spacing between grid points, and more points. Default is 30.
#' @param overlay logical, should an overlay contour be drawn? Default is TRUE
#' @param overlay_res the overlay resolution. Default is 60.
#' @param labels logical, should cluster labels be plotted? Default is TRUE.
#'
#' @returns A plot where lineages are represented as metro lines.
#'
#' @importFrom slingshot slingLineages
#' @importFrom ggplot2 .data aes geom_segment geom_point geom_polygon unit labs
#' @importFrom ggplot2 scale_color_manual theme_classic theme element_blank scale_size
#' @importFrom ggrepel geom_label_repel
#' @importFrom oveRlay makeOverlay
#'
#' @export

plotMetroMap <- function(sce,
                            coords_by,
                            dimred,
                            snap = TRUE,
                            res = 30,
                            overlay = TRUE,
                            overlay_res = 60,
                            include_cells = FALSE,
                            labels = FALSE) {

  lineages = slingLineages(sce)

  coords = cbind(tapply(reducedDim(sce, dimred)[,1], colData(sce)[,coords_by], median),
                 tapply(reducedDim(sce, dimred)[,2], colData(sce)[,coords_by], median))

  colnames(coords) = c("x", "y")
  coords = as.matrix(coords)

  nout = as.numeric(table(unlist(lineages)))

  lineages = lineages[lengths(lineages) > 1]

  linecoords = lapply(names(lineages), function(x) {
    from = lineages[[x]][seq_len(length(lineages[[x]]) - 1)]
    to = lineages[[x]][2 : length(lineages[[x]])]
    df = data.frame(from = from,
                    to = to,
                    x0 = coords[from, "x"],
                    x1 = coords[to, "x"],
                    y0 = coords[from, "y"],
                    y1 = coords[to, "y"])
    df$lineage = x
    return(df)
  })

  linecoords = do.call(rbind, linecoords)
  linecoords = linecoords[!duplicated(linecoords[,c("from", "to", "lineage")]),]


  # Colors
  colorpal = c("#346524", "#d04648", "#dad45e", "#d27d2c", "#30346d", "#140c1c",
               "#6dc2ca", "#442434", "#d2aa99", "#30346d", "#4e4a4e", "#757161",
               "#854c30", "#6daa2c")

  colorpal = colorpal[seq_along(lineages)]
  names(colorpal) = names(lineages)


  if(snap) {
    coordlist = snapDots(linecoords, res)
    linecoords = coordlist$segments
    grid = as.matrix(coordlist$grid)
    snc = snapToGrid(coords[,1:2], grid)
    rownames(snc) = rownames(coords)
    coords = snc
  }

  # Plot
  coords$nout = nout
  if(max(coords$nout) > 2) bks = seq(1, max(coords$nout), by = 2) else bks = c(0,1,2)

  if(!include_cells) {
  p = ggplot(linecoords, aes(col = .data[["lineage"]])) +
    geom_segment(aes(x = .data[["x0"]], xend = .data[["x1"]],
                     y = .data[["y0"]], yend = .data[["y1"]]),
                 size = 1.5) +
    geom_point(data = coords,
               mapping = aes(x = .data[["x"]],
                             y = .data[["y"]],
                             size = .data[["nout"]]),
               #size = 2,
               inherit.aes = FALSE, shape = 21, fill = "white") +
    scale_color_manual(values = colorpal) +
    scale_size(range = c(2, min(5, max(nout) + 2)), breaks = bks) +
    labs(size = "N. of lineages") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank()
          )
  } else {
    p = plot_UMAP(sce, umap_slot = dimred, color_by = coords_by)
    p = p +
      geom_segment(data = linecoords,
                   mapping = aes(x = .data[["x0"]], xend = .data[["x1"]],
                                 y = .data[["y0"]], yend = .data[["y1"]],
                                 col = .data[["lineage"]]),
                   inherit.aes = FALSE,
                   size = 1.5) +
      geom_point(data = coords,
                 mapping = aes(x = .data[["x"]], y = .data[["y"]],
                               size = .data[["nout"]]),
                 #size = 2,
                 inherit.aes = FALSE, shape = 21, fill = "white") +
      scale_color_manual(values = colorpal) +
      scale_size(range = c(2, min(5, max(nout) + 2)), breaks = bks) +
      labs(size = "N. of lineages")
  }

  if(overlay) {
    ov = makeOverlay(reducedDim(sce, dimred), res = overlay_res)
    p = p + geom_polygon(data = ov, aes(x = .data[["x"]],
                                        y = .data[["y"]],
                                        group = .data[["cluster_hole"]],
                                        subgroup = .data[["id_hole"]]),
                         size = 0.1,
                         col = "black",
                         fill = NA,
                         linetype = 2,
                         inherit.aes = FALSE)
  }

  if(labels) {
    clabs = as.data.frame(coords)
    clabs$label = rownames(clabs)
    p = p + geom_label_repel(data = clabs,
                             mapping = aes(x = .data[["x"]],
                                           y = .data[["y"]],
                                           label = .data[["label"]]),
                             inherit.aes = FALSE,
                             size = 3,
                             label.padding = unit(0.1, "lines"),
                             alpha = 0.7)
  }

  return(p)
}

#' @noRd
#' @importFrom pracma pdist2
snapToGrid <- function(dat, grid) {
  dmat = pdist2(X = dat, Y = grid)
  nn = unlist(apply(dmat, 1, which.min))
  if(any(duplicated(nn))) {
    reassign = which(duplicated(nn))
    for(i in seq_len(length(reassign))) {
      points = as.numeric(dmat[reassign[i],])
      points[which.min(points)] = max(points)*2
      new_point = which.min(points)
      while(new_point %in% nn) {
        points[which.min(points)] = max(points)*2
        new_point = which.min(points)
      }
      nn[reassign[i]] = new_point
    }
  }
  new_grid = as.data.frame(grid[nn,])
  return(new_grid)
}


#' @noRd
snapDots <- function(s, res) {

  sl = split(s, s$lineage)
  d = lapply(names(sl), function(x){
    p = makeLinesFromSegments(sl[[x]])
    p$lineage = x
    return(p)
  })

  names(d) = names(sl)

  stepsize = diff(range(s[,c("x0", "y0", "x1", "y1")]))/res
  grid = as.matrix(makeDotGrid(s[,c("x0", "y0", "x1", "y1")], stepsize))

  sn = lapply(names(d), function(x) {
    df = as.data.frame(snapToGrid(as.matrix(d[[x]][,c("x", "y")]), grid))
    df = makeSegmentsFromLine(df)
    df$lineage = x
    df = df[!duplicated(df),]
    return(df)
  })

  ss = do.call(rbind, sn)

  return(list("segments" = ss, "grid" = grid))
}

#' @noRd
makeSegmentsFromLine <- function(data) {

  df = list()
  for(i in seq_len(nrow(data)-1)) {
    x0 = data$x[i]
    x1 = data$x[i + 1]
    y0 = data$y[i]
    y1 = data$y[i + 1]
    df[[i]] = c(x0, x1, y0, y1)
  }
  df = as.data.frame(do.call(rbind, df))
  colnames(df) = c("x0", "x1", "y0", "y1")
  if("lineage" %in% colnames(data)) {
    df$lineage = unique(data$lineage)
  }
  return(df)
}

#' @noRd
makeLinesFromSegments <- function(segs) {

  xc = c(rbind(segs$x0, segs$x1))
  yc = c(rbind(segs$y0, segs$y1))

  df = data.frame("x" = xc, "y" = yc)
  df = df[!duplicated(df),]
  return(df)
}

#' @noRd
makeDotGrid = function(data, stepsize){
  maxstep = max(range(data)) + stepsize
  minstep = min(range(data)) - stepsize
  steps = seq(minstep, maxstep, by = stepsize)
  dots = as.data.frame(expand.grid(steps, steps))
  colnames(dots) = c("x", "y")
  return(dots)
}
