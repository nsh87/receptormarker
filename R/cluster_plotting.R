#' Plot of total within sum of squares for numbers of clusters
#' 
#' This function plots within sum of squares.
#'
#' @param clust_obj 
#' @param optimal 
#'
#' @return
#' @export
#'
#' @examples
wss_plot <- function(clust_obj, optimal = FALSE) {
  wss <- clust_obj['wss']
  start <- length(wss[is.na(wss)]) + 1 # We start the krange right after NA's
  end <- length(wss)
  krange <- start:end
  wss <- wss[!is.na(wss)]
  plot(krange, wss, type = "b", xlab = "# of Clusters", 
       ylab = "Within Sum of Squares", 
       main = "Within Sum of Squares by Cluster")
  if (optimal) {
    opti_clust <- clust_obj['k_best']
    points(opti_clust, wss[opti_clust], col = "red", pch = 1, cex = 3)
    legend("topright", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot of cluster gap statistic
#' 
#' This function plots gap analysis.
#'
#' @param clust_obj 
#' @param optimal 
#'
#' @return
#' @export
#'
#' @examples
clusGap_plot <- function(clust_obj, optimal = FALSE) {
  plot(clust_obj['clust_gap'], xlab = "# of Clusters", main = "Gap Analysis")
  if (optimal) {
    opti_clust <- clust_obj['k_best']
    gap_best <- clust_obj['clust_gap']['Tab'][[1]][opti_clust, 3]
    points(opti_clust, gap_best, col = "red", pch = 1, cex = 3)
    legend("topleft", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot of first two principal components.
#' 
#' This function plots the instances of \code{data} using the first two
#' principal components as the x and y axes, respectively.
#'
#' @param data 
#' @param clust_obj 
#' @param num_clust 
#'
#' @return
#' @export
#'
#' @examples
pca_plot <- function(data, clust_obj, num_clust) {
  pca <- princomp(data)
  clusters <- clust_obj['clust_model'][[num_clust]]['cluster']
  plot(pca$scores[,1:2], col = rainbow(num_clust)[clusters],
       xlab = "Principal Component 1",
       ylab = "Principal Component 2",
       main = "PCA Plot of Clusters")
}

#' Plot silhouette scores for a given clustering of data.
#' 
#' This function plots silhouette scores for each cluster in a given clustering
#' of data by leveraging the special way that the \code{\link{plot}} function
#' handles a \code{\link{silhouette}} object.
#'
#' @param clust_obj 
#' @param num_clust 
#'
#' @return
#' @export
#'
#' @examples
sil_plot <- function(clust_obj, num_clust) {
  sil <- clust_obj['sil'][[num_clust]]
  plot(sil)
}

#' Plot average silhouette widths for different numbers of clusters.
#' 
#' This function plots the average silhouette widths for all numbers of
#' clusters present in \code{clust_obj}.
#'
#' @param clust_obj 
#' @param optimal 
#'
#' @return
#' @export
#'
#' @examples
avg_sil_plot <- function(clust_obj, optimal = FALSE) {
  krange <- 1:length(clust_obj['sil_avg'])
  plot(krange, clust_obj['sil_avg'], type = "b", xlab = "# of Clusters",
       ylab = "Average Silhouette Width", 
       main = "Average Silhouette Width by Cluster")
  if (optimal) {
    opti_clust <- clust_obj['k_best']
    points(opti_clust, wss[opti_clust], col = "red", pch = 1, cex = 3)
    legend("bottomright", "Optimal Clusters", col = "red", pch = 1)
  }
}
