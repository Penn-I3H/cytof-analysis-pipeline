#' @useDynLib CL2, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider replace_na
#' @importFrom readr write_csv
#' @importFrom stringr str_split fixed str_replace str_remove
#' @importFrom tibble tibble as_tibble
#' @importFrom flowCore read.FCS flowFrame write.FCS
#' @importFrom flowDensity deGate
#' @importFrom KernSmooth bkde
#' @importFrom uwot umap
#' @importFrom RcppHNSW hnsw_knn
#' @importFrom ash bin2
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph graph_from_edgelist cluster_leiden V E E<- distances
#' @importFrom igraph membership components as_edgelist topological.sort bfs
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom sp point.in.polygon
#' @importFrom grDevices chull hcl hsv
#' @importFrom stats cor quantile sd setNames
#' @importFrom utils combn
NULL
