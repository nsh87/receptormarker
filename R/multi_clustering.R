#' Cluster with kmeans and find optimal k
#' 
#' This calls the function \code{\link[stats]{kmeans}} to perform kmeans 
#' clustering, but initializes multiple times. It chooses the best one for each
#' number of clusters using within sum of squares, and the best number of 
#' clusters overall based on the 26 methods used in the 
#' \code{\link[NbClust]{NbClust}} function.
#' 
#' @details The \code{multiClust} object is primarily intended for use with 
#'   the plotting functions in the \emph{Receptormarker} package (the
#'   \code{seealso} section contains a complete list). However, its elements can
#'   also be accessed directly for other uses.  
#'   All \code{\link[stats]{kmeans}} objects are available based on the
#'   \code{krange} passed into the function. That is, if one generated a
#'   \code{multiClust} object using a \code{krange} containing 3, and stored it
#'   in a variable called \code{clust_obj}, one could access the
#'   \code{\link[stats]{kmeans}} object for 3 clusters with the following call:  
#'   \code{clust_obj[["clust_model"]][[3]]}  
#'   
#'   This can be done similarly to access the \code{sil_avg}, \code{num_clust},
#'   \code{sil}, and \code{wss}. However, both \code{clust_gap} and
#'   \code{k_best} only contain single items. Thus, \code{clust_gap} can be
#'   accessed with the following call:  
#'   \code{clust_obj[["clust_gap"]]}  
#'   When selecting the values to pass in for each parameter, the defaults will
#'   likely perform well to start. However, a larger \code{krange} will generate
#'   more work for the function, so if certain values are definitely not going
#'   to be used and are not best, then consider not including them in
#'   \code{krange} for performance benefits. Reducing the number of \code{runs}
#'   will have performance benefits as well, though they may not be as
#'   pronounced. More runs increases the likelihood of finding the best
#'   clustering for a given number of clusters with a given algorithm.
#'   \code{iter.max} controls how many times the algorithm recalculates to try
#'   to converge on an ideal clustering. Without significant adjustments up or
#'   down in number, few changes should be observed here.  
#'   \strong{TIP:} When using the clustering portion of the
#'   \emph{Receptormarker} package, a reasonable workflow might be:  
#'   \enumerate{
#'     \item Create a \code{multiClust} object.
#'     \item Send this object to the \code{\link{wss_plot}},
#'       \code{\link{gap_plot}}, and \code{\link{avg_sil_plot}} with
#'       \code{optimal = TRUE} to see which number of clusters is recommended
#'       and compare it to these respective plots.
#'     \item Choose the number of clusters to be used.
#'     \item View the \code{\link{sil_plot}} using this chosen number of
#'       clusters to confirm the choice.
#'     \item View the \code{\link{pca_plot}} using the chosen number of clusters
#'       to see the shape of the them within the data based on using the first
#'       two principal components for axes.
#'     \item View the \code{\link{clust_boxplot}} using the chosen number of
#'       clusters to see the cluster membership of each feature (likely
#'       phenotypic markers). This should yield useful conclusions.
#'   } 
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
#'   \code{wss} A list of total within sum of squares for each of the
#'     \code{\link[stats]{kmeans}} objects stored in \code{clust_model}.
#'   \code{clust_gap} A \code{\link[cluster]{clusGap}} object that uses the
#'     highest number of \code{krange} for the \code{K.max} argument.
#'   \code{k_best} The optimal number of clusters based on 26 methods: "kl",
#'     "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
#'     "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2",
#'     "beale", "ratkowsky", "ball", "ptbiserial", frey", "mcclain", "dunn",
#'     "hubert", "sdindex", "dindex", "sdbw". Please see
#'     \code{\link[DbClust]{DbClust}} for further details on how each method
#'     works.
#' 
#' @export
#'
#' @seealso \code{\link{wss_plot}}, \code{\link{gap_plot}}, 
#'   \code{\link{pca_plot}}, \code{\link{sil_plot}}, \code{\link{avg_sil_plot}},
#'   \code{\link{clust_boxplot}}
#'
#' @examples
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
multi_clust <- function(d, krange = 2:15, iter.max = 300, runs = 10, 
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
  nb_best <- NbClust::NbClust(d, min.nc = krange[1],
                              max.nc = krange[length(krange)],
                              method = "average")
  best <- aggregate(nb_best$Best.nc[1,], by=list(nb_best$Best.nc[1,]), length)
  index <- which.max(best[[2]])
  km[["k_best"]] <- best[index, 1] 
  structure(km, class = "multiClust")
}

#' @title Validate and sort an argument that is a range of integers
#' @description An internal function that raises an error if the argument is not
#' a positive, non-duplicate range of integers beginning > 1. It is sorted into
#' ascending order.
#' @param n_range An item to be checked to make sure it is a valid range.
#' @keywords internal
validate_sort_range <- function(n_range) {
  if (class(n_range) != "integer" || length(n_range) <= 1) {
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
