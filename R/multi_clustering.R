#' Cluster with kmeans and find optimal k
#' 
#' This calls the function \code{\link[stats]{kmeans}} to perform kmeans 
#' clustering, but initializes multiple times. It chooses the best one for each
#' number of clusters using within sum of squares, and the best number of 
#' clusters overall based on average silhouette width from the 
#' \code{\link[cluster]{silhouette}} function.
#'
#' @param d A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns).
#' @param krange An integer vector. Numbers of clusters which are to be compared
#'   by the average silhouette width criterion. Note: \code{krange} should not
#'   include 1 since silhouette scores are not defined there.
#' @param iter.max An integer. The maximum number of iterations used by kmeans
#'   to find its centers.
#' @param runs An integer. The number of starts of the k-means algorithm.
#' @param method A string. Either "kmeans", "hclust", or "both". Currently, it
#'   does nothing, but functionality will be added in the future.
#' @param ... Further arguments to be passed to \code{\link[stats]{kmeans}}.
#'
#' @return \code{multi_clust} returns an object of class "multiClust"
#'   that can be used for multiple plots. It is a list with at least the
#'   following components:
#'   \code{clust_model} A list of \code{\link[stats]{kmeans}} objects for each
#'     number of cluster centers requested.
#'   \code{sil_avg} A list of average silhouette scores for each number of
#'     cluster centers requested.
#'   \code{num_clust} A list of the each number of clusters used.
#'   \code{sil} A list of \code{\link[cluster]{silhouette}} objects for each
#'     number of clusters.
#'   \code{clust_gap} A \code{\link[cluster]{clusGap}} object that uses the
#'     highest number of \code{krange} for the \code{K.max} argument.
#'   \code{k_best} The optimal number of clusters based on silhouette score.
#' 
#' @export
#'
#' @seealso \code{\link[stats]{kmeans}}, \code{\link[cluster]{silhouette}}, 
#'   \code{\link[cluster]{clusGap}}
#'
#' @examples
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
multi_clust <- function(d, krange = 2:10, iter.max = 300, runs = 10, 
                            method = "kmeans", ...) {
  validate_not_null(list(d = d, krange = krange, iter.max = iter.max, 
                         runs = runs, method = method))
  validate_num_data(d)
  krange <- validate_sort_range(krange) # returns sorted krange by ascending
  validate_pos_num(list(iter.max = iter.max, runs = runs))
  d_dist <- dist(d)
  km <- list(clust_model = NULL, sil_avg = NULL, num_clust = NULL, sil = NULL,
             clust_gap = NULL, wss = NULL, k_best = NULL)
  for (k in krange) {
    min_wss <- Inf
    km_opt <- NULL
    for (i in 1:runs) {
      kmm <- stats::kmeans(d, k, iter.max = iter.max, nstart = 10)
      swss <- kmm[["tot.withinss"]]
      if (swss < min_wss) {
        min_wss <- swss
        km_opt <- kmm
      }
    }
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
  km[["k_best"]] <- which.max(km[["sil_avg"]])
  structure(km, class = "multiClust")
}

#' @title Validate and sort an argument that is a range of integers
#' @description An internal function that raises an error if the argument is not
#' a positive, non-duplicate range of integers beginning > 1. It is sorted into
#' ascending order.
#' @param n_range An item to be checked to make sure it is a valid range.
#' @keywords internal
validate_sort_range <- function(n_range) {
  if (class(n_range) != "integer" || length(n_range <= 1)) {
    stop("The argument 'krange' must be a range of integers.", call.=FALSE)
  } else if (any(n_range <= 1)) {
    stop("The argument 'krange' must contain only values greater than one.", 
         call.=FALSE)
  } else if (anyDuplicated(n_range) != 0) {
    stop("The argument 'krange' must not contain any duplicate values.",
         call.=FALSE)
  }
  sort(n_range)
}

#' @title Validate that an argument contains positive integers
#' @description An internal function that raises an error if the argument does
#' not contain positive integers.
#' @param n A named list of items to be checked.
#' @keywords internal
validate_pos_num <- function(n) {
  lapply(1:length(n), 
         function(i) {
           if (class(n[[i]]) != "integer" || length(n[[i]]) > 1) {
             err <- paste0("The argument '", names(n)[[i]],
                           "' must be a positive integer", collapse="")
             stop(err, call.=FALSE)
           }
         })

#' @title Validate that an object is of class \emph{multiClust}
#' @description An internal function that raises an error if the argument is not
#' of class \emph{multiClust}.
#' @param clust_obj An item to be checked for class membership.
#' @keywords internal
validate_multi_clust <- function(clust_obj) {
  if (class(clust_obj) != "multiClust") {
    stop("The argument 'clust_obj' must be an object of class 'multiClust",
         call.=FALSE)
  }
}
