#' Title
#'
#' @param data 
#' @param krange 
#' @param iter.max 
#' @param runs 
#' @param ... 
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
  structure(km)
}
