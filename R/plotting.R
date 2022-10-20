


#' Plot cluster modularity
#'
#' Plots a heatmap showing pairwise cluster modularity from a SCE object
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "modularity_SNN_100"
#' @param type character, one of "heatmap" or "graph" for different plot types
#'
#' @return a heatmap of pairwise modularity
#'
#' @importFrom pheatmap pheatmap
#' @importFrom colorspace sequential_hcl
#' @importFrom igraph graph_from_adjacency_matrix E layout_with_lgl
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export


plotModularity <- function(sce, name, type = "heatmap") {
  name = paste0("modularity_", name)
  modmat = log2(metadata(sce)[[name]] + 1)
  if(type == "heatmap"){
  pheatmap(mat = modmat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = sequential_hcl(palette = "Sunset", n = 25),
           main = "log2(modularity ratio + 1)",
           border_color = NA
  )
  } else if(type == "graph") {

    cgr <- graph_from_adjacency_matrix(modmat,
                                       mode="upper",
                                       weighted=TRUE,
                                       diag=FALSE)

    set.seed(420)
    plot(cgr,
         edge.width=E(cgr)$weight*5,
         layout=igraph::layout_with_lgl,
         main = name)
  }
}


#' Plot approximate silhouette widths
#'
#' Plots a beesswarm plot showing the approximate width of each cell within
#' each cluster
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "silhouette_SNN_100"
#'
#' @return a beeswarm plot of silhouette widhts
#'
#' @importFrom ggplot2 ggplot aes_string theme_bw
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export

plotSilhouette <- function(sce, name) {
  name = paste0("silhouette_", name)
  ggplot(metadata(sce)[[name]], aes_string(x="cluster", y="width", colour="closest")) +
    geom_quasirandom(method="smiley") +
    theme_bw()
}
