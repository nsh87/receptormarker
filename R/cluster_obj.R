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
  
  km <- lapply(krange, 
               function(k) {
                 sil_max <- 0
                 km_opt <- NULL
                 for (i in 1:runs) {
                   kmm <- kmeans(data, k, iter.max = iter.max, ...)
                   sil_sum <- summary(silhouette(kmm$cluster, data_dist))
                   sil_avg <- sil_sum$avg.width
                   if (sil_avg > sil_max) {
                     sil_max <- sil_avg
                     km_opt <- kmm
                   }
                 }
                 list(km_obj = km_opt, sil_avg = sil_max)
               })
  k_best <- which.max(km$sil_avg)
}
