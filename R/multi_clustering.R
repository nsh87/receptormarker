#' @title Cluster using kmeans over multiple k and find optimal k
#' @description This calls the function \code{\link[stats]{kmeans}} to perform
#' kmeans clustering, but initializes multiple times. It chooses the best
#' clustering for each number of clusters (there are multiple runs for each
#' k) based on within sum of squares. The best overall k (i.e. best number of
#' clusters to use for the data set) is then chosen by running the data set
#' through a number of different cluster estimation algorithms. Majority rule
#' of all the algorithms determines the optimal k. For example, if 12 algorithms
#' suggest that k=9 and 7 algorithms suggest that k=6 then the optimal k will
#' be 9.
#' @details The \code{\link{multiClust-class}} object returned is primarily
#' intended for use with the plotting functions in the \emph{receptormarker}
#' package (the See Also section contains a complete list). However, its
#' elements can also be accessed directly for other uses. See the
#' \code{\link{multiClust-class}} documentation for details.
#'
#' When choosing the values to pass in for each parameter, the defaults will
#' likely perform well to start. However, a large \code{krange} will generate
#' more work for the function, so if certain values are definitely not going
#' to be used and are not best, then consider not including them in
#' \code{krange} for performance benefits. Reducing the number of \code{runs}
#' will have performance benefits as well, though they may not be as
#' pronounced. More runs increases the likelihood of finding the best
#' clustering for a given number of clusters, k, with a given algorithm.
#' \code{iter.max} controls how many times the algorithm recalculates to try
#' to converge on an ideal clustering. Without significant adjustments up or
#' down in number, few changes in performance should be observed here.  
#' 
#' \strong{TIP:} When using the clustering portion of the
#' \emph{receptormarker} package, a reasonable workflow might be:  
#' \enumerate{
#'   \item Create a \code{\link{multiClust-class}} object using this function.
#'   \item Send that object to the functions \code{\link{wss_plot}},
#'     \code{\link{gap_plot}}, and \code{\link{avg_sil_plot}} with
#'     \code{optimal = TRUE} to see which number of clusters is recommended
#'     and compare it to these respective plots.
#'   \item Choose the number of clusters to be used.
#'   \item View the \code{\link{sil_plot}} using your chosen number of
#'     clusters to confirm the choice.
#'   \item View the \code{\link{pca_plot}} using your chosen number or range
#'     of clusters to see the shape of the data using its first two principal
#'     components.
#'   \item View the \code{\link{clust_boxplot}} using your chosen number of
#'     clusters to see the cluster membership of each feature (likely
#'     phenotypic markers). This should yield useful conclusions.
#' } 
#' @param d A numeric matrix of data, or an object that can be coerced to
#' such a matrix (such as a numeric vector or a data frame with all numeric
#' columns).
#' @param krange An integer vector describing the range of k (numbers of
#' clusters) which are to be compared. Note: \code{krange} should not
#' include 1 since silhouette scores and clusters are not defined at this value.
#' @param iter.max An integer representing the maximum number of iterations to
#' be find cluster centers when running kmeans. Only used when the \code{method}
#' parameter is set to "kmeans".
#' @param runs An integer representing the number of starts of the k-means
#' algorithm for each k. Only used when the \code{method} parameter is set to
#' "kmeans".
#' @param method Either "kmeans" or "exhaustive". If "kmeans", then the average
#' silhouette score will be used to estimate the optimal k using KMeans
#' clustering. If "exhausting", then \code{\link{NbClust}} will be used to
#' estimate optimal k, although this takes significantly longer and is more
#' error prone since the data is run through upwards of 20 algorithms for
#' clustering; consider this method to be in beta.
#' @param index The index to be calculated for determining the optimal k, passed
#' to \code{\link{NbClust}} (see its documentation for details on this
#' parameter). Only used when \code{method} is set to "exhaustive".
#' @param ... Further arguments to be passed to \code{\link[stats]{kmeans}}.
#' @return An object of class \code{\link{multiClust-class}}. This object can
#' be used to create several plots (refer to the See Also section) that aid
#' in determining the optimal k for the data set. It also contains statistics
#' about each clustering in \code{krange}.
#' @export
#' @seealso \code{\link{wss_plot}}, \code{\link{gap_plot}}, 
#' \code{\link{pca_plot}}, \code{\link{sil_plot}}, \code{\link{avg_sil_plot}},
#' \code{\link{clust_boxplot}}, \code{\link{multiClust-class}}
#' @examples
#' # Common data set
#' library(datasets)
#' iris_clust <- multi_clust(iris[, 1:4])
#' # Domain data set
#' data(fluidigm)
#' fluidigm_clust <- multi_clust(fluidigm[1:40, ])
multi_clust <- function(d, krange = 2:15, iter.max = 500, runs = 10, 
                        method = "kmeans", index = "alllong", ...) {
  validate_not_null(list(d = d, krange = krange, iter.max = iter.max, 
                         runs = runs, method = method, index = index))
  validate_k_range(krange)
  validate_pos_num(list(iter.max = iter.max, runs = runs))
  validate_num_data(d)
  validate_multi_clust_method(method)
  krange <- sort(krange)
  bool <- is_boolean(d)
  d <- d[complete.cases(d), ]
  d_dist <- dist(d)
  km <- list(clust_model = NULL, sil_avg = NULL, num_clust = NULL, sil = NULL,
             clust_gap = NULL, wss = NULL, k_best = NULL)
  # Perform KMeans clustering for every k in krange
  for (k in krange) {
    km_opt <- stats::kmeans(d, k, iter.max = iter.max, nstart = runs, ...)
    min_wss <- km_opt[["tot.withinss"]]
    sil <- cluster::silhouette(km_opt[["cluster"]], d_dist)
    sil_sum <- summary(sil)
    sil_avg <- sil_sum[["avg.width"]]
    
    km[["clust_model"]][[k]] <- km_opt
    km[["sil_avg"]][[k]] <- sil_avg
    km[["num_clust"]][[k]] <- k
    km[["sil"]][[k]] <- sil
    km[["wss"]][[k]] <- min_wss
  }
  km[["clust_gap"]] <- cluster::clusGap(d, kmeans, 
                                        K.max = length(km[["clust_model"]]), 
                                        B = 15, verbose = FALSE)
  # Estimate K using the average silhouette score if method == 'kmeans'
  if (method == "kmeans") {
    km[["k_best"]] <- which.max(km[["sil_avg"]])
  } else {
  # Estimate K using NbClust if method == 'exhaustive'
    if (bool) {
      distance <- "binary"
    } else {
      distance <- "euclidean"
    }
    tryCatch({
      nb_best <- suppressWarnings(suppressMessages(
                NbClust(d,
                min.nc = krange[1],
                index = index,
                max.nc = krange[length(krange)],
                distance = distance,
                method = "average")))
    },
    error = function(e) {
      if (grepl("computationally singular", e)) {
        stop("There are not enough rows of data to evaluate for clustering.",
             call. = FALSE)
      } else {
        stop(e, call. = FALSE)
      }
    }
    )
    best <- aggregate(nb_best[["Best.nc"]][1, ], 
                      by = list(nb_best[["Best.nc"]][1, ]), length)
    idx <- which.max(best[[2]])
    km[["k_best"]] <- best[idx, 1]
  }
  new("multiClust", clust_model=km[["clust_model"]], sil_avg=km[["sil_avg"]],
      num_clust=km[["num_clust"]], sil=km[["sil"]], clust_gap=km[["clust_gap"]],
      wss=km[["wss"]], k_best=km[["k_best"]])
}
