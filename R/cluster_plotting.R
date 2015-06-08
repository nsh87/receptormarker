#' Plot of total within sum of squares for numbers of clusters
#' 
#' This function plots.
#'
#' @param clust_obj 
#' @param optimal 
#'
#' @return
#' @export
#'
#' @examples
wss_plot <- function(clust_obj, optimal = FALSE) {
  krange <- 1:length(clust_obj$clust_model)
  wss <- vapply(krange, 
                function(k) {
                  if (is.na(clust_obj$clust_model[[k]])) {
                    wss <- 0
                  } else {
                    wss <- sum(clust_obj$clust_model[[k]]$withinss)
                  }
                  wss
                },
                numeric(1))
  plot(krange, wss,type="b", xlab="# of Clusters", ylab="Within Sum of Squares",
       main="Within Sum of Squares by Cluster")
  if (optimal) {
    opti_clust <- clust_obj$k_best
    points(opti_clust, wss[opti_clust], col = "red", pch = 1, cex = 3)
    legend("topright", "Optimal Clusters", col = "red", pch = 1)
  }
}
