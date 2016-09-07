#' Calculate total within sum of squares for a clustering and distance object
#'
#' This is an internal function that returns the total sum of squares for a
#' given data set by passing it the distance object and cluster assignments. It
#' should only be used when the clustering object does not contain this
#' information (e.g. an \code{hclust} object).
#' @param dist_obj The distance object for the data set that has been
#'   clustered.
#' @param clusters The cluster assignments the distance object has. This should
#'   be a vector.
#'
#' @return An integer value indicating the total within sum of squares
#' @keywords internal
totwss <- function(dist_obj, clusters){
  totwss <- 0
  for (cluster in unique(clusters)) {
    indices <- which(cluster == clusters)
    n <- length(indices)
    if (n > 1) {
      clust_dist_obj <- as.dist(as.matrix(dist_obj)[indices, indices])
      wss <- mean(clust_dist_obj ^ 2) / 2 * (n - 1)
    } else if (n == 1) {
      wss <- 0
    } else {
      wss <- NA
    }
    totwss <- totwss + wss
  }
  return (totwss)
}
