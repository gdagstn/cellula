


#' Plot cluster modularity
#'
#' Plots a heatmap showing pairwise cluster modularity from a SCE object
#'
#' @param sce a SingleCellExperiment object
#' @param name character, the name in the metadata slot of the SCE object, e.g.
#'     "modularity_SNN_100"
#'
#' @return a heatmap of pairwise modularity
#'
#' @importFrom pheatmap pheatmap
#' @importFrom colorspace sequential_hcl
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export


plotModularity <- function(sce, name) {
  name = paste0("modularity_", name)
  pheatmap(mat = log2(metadata(sce)[[name]] + 1),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = sequential_hcl(palette = "Sunset", n = 25),
           main = "log2(modularity ratio + 1)",
           border_color = NA
  )
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
