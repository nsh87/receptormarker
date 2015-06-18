#' Plot of total within sum of squares for numbers of clusters
#' 
#' This function plots within sum of squares. It provides a wrapper on 
#' \code{\link[graphics]{plot}} by extracting \code{wss} from \code{clust_obj}
#' and plotting it.
#'
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{wss}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#' 
#' @seealso \code{\link[graphics]{plot}}, \code{\link{cluster_optimal}}
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with wss_plot
#' wss_plot(iris_cluster, optimal = TRUE)
wss_plot <- function(clust_obj, optimal = FALSE, ...) {
  wss <- clust_obj[['wss']]
  start <- length(wss[is.na(wss)]) + 1 # We start the krange right after NA's
  end <- length(wss)
  krange <- start:end
  wss <- wss[!is.na(wss)]
  plot(krange, wss, type = "b", xlab = "# of Clusters", 
       ylab = "Within Sum of Squares", 
       main = "Within Sum of Squares by Cluster",
       ...)
  if (optimal) {
    opti_clust <- clust_obj[['k_best']]
    index <- opti_clust - (start - 1) # Account for removal of NA's
    points(opti_clust, wss[index], col = "red", pch = 1, cex = 3)
    legend("topright", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot of cluster gap statistic
#' 
#' This function plots gap analysis. It provides a wrapper on 
#' \code{\link[graphics]{plot}} by extracting \code{clust_gap}, a 
#' \code{\link[cluster]{clusGap}} object, and leveraging the custom plotting
#' method that it has.
#'
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{clust_gap}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link[graphics]{plot}}, \code{\link{cluster_optimal}},
#'   \code{\link[cluster]{clusGap}}
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with clusGap_plot
#' clusGap_plot(iris_cluster, optimal = TRUE)
clusGap_plot <- function(clust_obj, optimal = FALSE, ...) {
  plot(clust_obj[['clust_gap']], xlab = "# of Clusters", main = "Gap Analysis",
       ...)
  if (optimal) {
    opti_clust <- clust_obj[['k_best']]
    gap_best <- clust_obj[['clust_gap']]['Tab'][[1]][opti_clust, 3]
    points(opti_clust, gap_best, col = "red", pch = 1, cex = 3)
    legend("topleft", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot of first two principal components.
#' 
#' This function plots the instances of \code{data} using the first two
#' principal components as the x and y axes, respectively. It represents each
#' cluster with a different color so that one can understand the distribution of
#' the clusters based on the first two components. It is a wrapper on the 
#' \code{\link[graphics]{plot}} function.
#'
#' @param data A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns). Note: This should be the same one used to generate 
#'   \code{clust_obj}.
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{\link{cluster_optimal}} object.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link[graphics]{plot}}, \code{\link{cluster_optimal}},
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with pca_plot
#' pca_plot(iris[, 1:4], iris_cluster, num_clust = 3)
pca_plot <- function(data, clust_obj, num_clust, ...) {
  pca <- princomp(data)
  clusters <- clust_obj[['clust_model']][[num_clust]][['cluster']]
  plot(pca[['scores']][,1:2], col = rainbow(num_clust)[clusters],
       xlab = "Principal Component 1",
       ylab = "Principal Component 2",
       main = "PCA Plot of Clusters",
       ...)
}

#' Plot silhouette scores for a given clustering of data.
#' 
#' This function plots silhouette scores for each cluster in a given clustering
#' from \code{clust_obj} by extracting a given \code{sil} from it and leveraging
#' the special way that the \code{\link[graphics]{plot}} function handles a
#' \code{\link[cluster]{silhouette}} object.
#'
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{\link{cluster_optimal}} object.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link[graphics]{plot}}, \code{\link{cluster_optimal}},
#' \code{\link[cluster]{silhouette}}
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with sil_plot
#' sil_plot(iris_cluster, num_clust = 3)
sil_plot <- function(clust_obj, num_clust, ...) {
  sil <- clust_obj[['sil']][[num_clust]]
  plot(sil, main = "Silhouette Plot", ...)
}

#' Plot average silhouette widths for different numbers of clusters.
#' 
#' This function plots the average silhouette widths for all numbers of
#' clusters present in \code{clust_obj}. It provides a wrapper on 
#' \code{\link[graphics]{plot}} by extracting \code{sil_avg} from 
#' \code{clust_obj} and plotting it.
#'
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{clust_gap}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link{plot}}
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link[graphics]{plot}}, \code{\link{cluster_optimal}},
#'   \code{\link[cluster]{silhouette}}
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with avg_sil_plot
#' avg_sil_plot(iris_cluster, optimal = TRUE)
avg_sil_plot <- function(clust_obj, optimal = FALSE, ...) {
  sil <- clust_obj[['sil_avg']]
  start <- length(sil[is.na(sil)]) + 1 # We start the krange right after NA's
  end <- length(sil)
  krange <- start:end
  sil <- sil[!is.na(sil)]
  plot(krange, sil, type = "b", xlab = "# of Clusters",
       ylab = "Average Silhouette Width", 
       main = "Average Silhouette Width by Cluster",
       ...)
  if (optimal) {
    opti_clust <- clust_obj[['k_best']]
    index <- opti_clust - (start - 1) # Account for removal of NA's
    points(opti_clust, sil[index], col = "red", pch = 1, cex = 3)
    legend("bottomright", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot cluster membership for each feature.
#' 
#' This function plots cluster membership for each feature of \code{data} using
#' box plots. It utlizes \code{facet_wrap} in the \code{\link[ggplot2]{qplot}}
#' function and the desired \code{clust_model} from the argument 
#' \code{clust_obj} to show the box plots of each feature after using the
#' \code{\link[reshape2]{melt}} function to get the correct form of \code{data}.
#'
#' @param data A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns). Note: This should be the same one used to generate 
#'   \code{clust_obj}.
#' @param clust_obj A \code{\link{cluster_optimal}} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{\link{cluster_optimal}} object.
#' @param ... Further arguments to be passed to the \code{\link{qplot}} 
#'   function (besides \code{xlab}, \code{ylab}).
#'
#' @export
#'
#' @seealso \code{\link[ggplot2]{qplot}}, \code{\link{cluster_optimal}},
#'   \code{\link[reshape2]{melt}}
#'
#' @examples
#' # First, create a cluster_optimal object
#' library(datasets)
#' iris_cluster <- cluster_optimal(iris[, 1:4])
#' # Second, use object with clust_boxplot
#' clust_boxplot(iris[, 1:4], iris_cluster, num_clust = 3)
clust_boxplot <- function(data, clust_obj, num_clust, ...) {
  meas_vars <- colnames(data)
  data['cluster'] <- clust_obj[['clust_model']][[num_clust]][['cluster']]
  m <- reshape2::melt(data, id.vars = "cluster", measure.vars = meas_vars)
  ggplot2::qplot(x = as.factor(cluster), y = value, data = m, geom = "boxplot", 
                 fill = as.factor(cluster), xlab = NULL, 
                 ylab = "Relative expression level", ...) + 
    facet_wrap(~variable) + scale_fill_discrete(name = "Cluster")
}
