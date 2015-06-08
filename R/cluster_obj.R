#' Cluster with kmeans and find optimal k
#' 
#' This calls the function \code{\link{kmeans}} to perform kmeans clustering,
#' but initializes multiple times. It chooses the best one for each number of
#' clusters, and the best number of clusters overall, based on average 
#' silhouette width using the \code{\link{silhouette}} function.
#'
#' @param data A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns).
#' @param krange An integer vector. Numbers of clusters which are to be compared
#'   by the average silhouette width criterion. Note: \code{krange} should not
#'   include 1 since silhouette scores are not defined there.
#' @param iter.max An integer. The maximum number of iterations allowed.
#' @param runs An integer. Number of starts of the k-means algorithm.
#' @param ... further arguments to be passed to \code{\link{kmeans}}.
#'
#' @return \code{cluster_obj} returns an object of class "cluster_obj" that can
#'   be used for multiple plots. It is a list with at least the following
#'   components:
#'   \code{clust_model} A list of \code{\link{kmeans}} objects for each number
#'     of cluster centers requested.
#'   \code{sil_avg} A list of average silhouette scores for each number of
#'     cluster centers requested.
#'   \code{num_clust} A list of the each number of clusters used.
#'   \code{sil} A list of \code{\link{silhouette}} objects for each number of
#'     clusters.
#'   \code{clust_gap} A \code{\link{clusGap}} object that uses the highest
#'     number of \code{krange} for the \code{K.max} argument.
#'   \code{k_best} The optimal number of clusters based on silhouette score.
#' 
#' @import 
#' 
#' @export
#'
#' @seealso \code{\link{kmeans}}, \code{\link{silhouette}}
#'
#' @examples
#' library(datasets)
#' iris_cluster <- cluster_obj(iris[, 1:4])
cluster_obj <- function(data, krange = 2:10, iter.max = 100, runs=100, ...) {
  if (1 %in% krange) stop("The entire range for # of clusters is to be > 1.")
  data_dist <- dist(data)
  km <- list(clust_model = NULL, sil_avg = NULL, num_clust = NULL, sil = NULL,
             clust_gap = NULL, k_best = NULL)
  for (k in krange) {
    sil_max <- 0
    sil_avg_max <- 0
    km_opt <- NULL
    for (i in 1:runs) {
      kmm <- kmeans(data, k, iter.max = iter.max)
      sil <- cluster::silhouette(kmm$cluster, data_dist)
      sil_sum <- summary(sil)
      sil_avg <- sil_sum$avg.width
      if (sil_avg > sil_avg_max) {
        sil_avg_max <- sil_avg
        km_opt <- kmm
        sil_max <- sil
      }
    }
    km$clust_model[[k]] <- km_opt
    km$sil_avg[[k]] <- sil_avg_max
    km$num_clust[[k]] <- k
    km$sil[[k]] <- sil_max
  }
  km$clust_gap <- cluster::clusGap(data, kmeans, 
                                   K.max = length(km$clust_model), B = 15,
                                   verbose = FALSE)
  km$k_best <- which.max(km$sil_avg)
  structure(km, class = "cluster_obj")
}
