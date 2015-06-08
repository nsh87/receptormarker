#' Cluster with kmeans and find optimal k
#' 
#' This calls the function \code{\link{kmeans}} to perform kmeans clustering,
#' but initializes multiple times and chooses the best one for each number of
#' clusters based on average silhouette width using the \code{\link{silhouette}}
#' 
#'
#' @param data A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns).
#' @param krange An integer vector. Numbers of clusters which are to be compared
#'   by the average silhouette width criterion.
#' @param iter.max An integer. The maximum number of iterations allowed.
#' @param runs An integer. Number of starts of the k-means algorithm.
#' @param ... further arguments to be passed to \code{\link{kmeans}}.
#'
#' @return 
#' 
#' @import cluster
#' 
#' @export
#'
#' @examples
cluster_obj <- function(data, krange = 2:10, iter.max = 100, runs=100, ...) {
  if (1 %in% krange) stop("The entire range for # of clusters is to be > 1.")
  data_dist <- dist(data)
  km <- list(km_obj = NULL, sil_avg = NULL, num_clusters = NULL, sil = NULL)
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
    km$km_obj[[k]] <- km_opt
    km$sil_avg[[k]] <- sil_avg_max
    km$num_clusters[[k]] <- k
    km$sil[[k]] <- sil_max
  }
  km$k_best <- which.max(km$sil_avg)
  structure(km, class = "cluster_obj")
}
