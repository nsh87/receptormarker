#' Plot of total within sum of squares for numbers of clusters
#' 
#' This function plots within sum of squares. It provides a wrapper on 
#' \code{\link[graphics]{plot}} by extracting \code{wss} from \code{clust_obj}
#' and plotting it.
#' 
#' @details This function is intended to be used as part of deciding the ideal
#'   number of clusters to use for analysis. Within sum of squares (WSS) tells
#'   one how tightly the data points within a cluster are clustered. Typically,
#'   one looks for an "elbow" at which the WSS drops significantly. Sometimes,
#'   the elbow in the WSS still does not contain enough clusters to understand
#'   the data well. This is where other visualizations can be useful.  
#'   Choosing \code{optimal = TRUE} will circle the optimal number of clusters
#'   based on average silhouette width. See the \emph{Details} section for the
#'   \code{\link{multi_clust}} function and view the \emph{TIP} for a suggested
#'   workflow.
#'
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{wss}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#' 
#' @seealso \code{\link{multi_clust}}, \code{\link{gap_plot}}, 
#'   \code{\link{pca_plot}}, \code{\link{sil_plot}}, \code{\link{avg_sil_plot}},
#'   \code{\link{clust_boxplot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with wss_plot
#' wss_plot(iris_cluster, optimal = TRUE)
wss_plot <- function(clust_obj, optimal = FALSE, ...) {
  validate_not_null(list(clust_obj = clust_obj, optimal = optimal))
  validate_true_false(list(optimal = optimal))
  validate_multi_clust(clust_obj)
  wss <- clust_obj[["wss"]]
  start <- length(wss[is.na(wss)]) + 1 # We start the krange right after NA's
  end <- length(wss)
  krange <- start:end
  wss <- wss[!is.na(wss)]
  plot(krange, wss, type = "b", xlab = "# of K Clusters",
       ylab = "Within Sum of Squares", 
       main = "Within Sum of Squares by K Clusters",
       ...)
  if (optimal) {
    opti_clust <- clust_obj[["k_best"]]
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
#' @details This function is intended to be used as part of deciding the ideal
#'   number of clusters to use for analysis. The cluster gap statistic tells
#'   one the goodness of a given clustering using the gap statistic. Typically,
#'   one looks for the peak before which the gap rises significantly. See
#'   \code{\link[cluster]{clusGap}} for further details on its calculation.
#'   Sometimes, the peak in the gap still does not contain enough clusters to
#'   understand the data well. This is where other visualizations can be useful.  
#'   Choosing \code{optimal = TRUE} will circle the optimal number of clusters
#'   based on average silhouette width. See the \emph{Details} section for the
#'   \code{\link{multi_clust}} function and view the \emph{TIP} for a suggested
#'   workflow.
#'
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{clust_gap}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link{multi_clust}}, \code{\link{wss_plot}}, 
#'   \code{\link{pca_plot}}, \code{\link{sil_plot}}, \code{\link{avg_sil_plot}},
#'   \code{\link{clust_boxplot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with gap_plot
#' gap_plot(iris_cluster, optimal = TRUE)
gap_plot <- function(clust_obj, optimal = FALSE, ...) {
  validate_not_null(list(clust_obj = clust_obj, optimal = optimal))
  validate_true_false(list(optimal = optimal))
  validate_multi_clust(clust_obj)
  plot(clust_obj[["clust_gap"]], xlab = "# of K Clusters",
       main = "Gap Analysis by K Clusters", ...)
  if (optimal) {
    opti_clust <- clust_obj[["k_best"]]
    gap_best <- clust_obj[["clust_gap"]]["Tab"][[1]][opti_clust, 3]
    points(opti_clust, gap_best, col = "red", pch = 1, cex = 3)
    legend("topleft", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Plot of first two principal components.
#' 
#' This function plots the instances of \code{data} using the first two
#' principal components as the x and y axes, respectively. These components are
#' calculated using \code{\link[stats]{princomp}}. It represents each cluster
#' with a different color so that one can understand the distribution of the
#' clusters based on the first two components. It is a wrapper on the
#' \code{\link[graphics]{plot}} function.
#' 
#' @details This function is intended to be used to view the shape of a given
#'   clustering to use for analysis. Principal component analysis looks for the
#'   way to represent the data such that the most variance is explained using a
#'   linear combination of the features. See \code{\link[stats]{princomp}} for
#'   further details on its calculation. Here it is used strictly for its
#'   ability to represent the data relatively well in just a couple dimensions.  
#'   One may pass in the \code{k_best} element from the \code{multiClust} object
#'   for \code{num_clust}, or use a different value. This is where other
#'   visualizations can be useful.  
#'   See the \emph{Details} section for the \code{\link{multi_clust}} function
#'   and view the \emph{TIP} for a suggested workflow.
#'
#' @param d A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns). Note: This should be the same one used to generate 
#'   \code{clust_obj}.
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{multiClust} object.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link{multi_clust}}, \code{\link{wss_plot}}, 
#'   \code{\link{gap_plot}}, \code{\link{sil_plot}},
#'   \code{\link{avg_sil_plot}}, \code{\link{clust_boxplot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with pca_plot
#' pca_plot(iris[, 1:4], iris_cluster, num_clust = 3)
pca_plot <- function(d, clust_obj, num_clust, ...) {
  validate_num_data(d)
  validate_multi_clust(clust_obj)
  validate_pos_num(list(num_clust = num_clust))
  pca <- stats::princomp(d)
  sdev <- pca[["sdev"]]
  prop_var <- sdev ^ 2 / sum(sdev ^ 2)
  main <- paste0("PCA Plot (", round(sum(prop_var[1:2]) * 100), "% Variance)")
  clusters <- clust_obj[["clust_model"]][[num_clust]][["cluster"]]
  clust_colors <- rainbow(num_clust)[clusters]
  par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
  plot(pca[["scores"]][, 1:2], col = clust_colors,
       xlab = "Principal Component 1",
       ylab = "Principal Component 2",
       main = main,
       ...)
  legend("topright",
         inset = c(-0.3, 0),
         legend = 1:num_clust,
         lty = rep(1, num_clust),
         col = clust_colors,
         cex = 0.6,
         title = "Cluster")
  par(mar = c(5.1, 4.1, 4.1, 4.1))
}

#' Plot silhouette scores for a given clustering of data.
#' 
#' This function plots silhouette scores for each cluster in a given clustering
#' from \code{clust_obj} by extracting a given \code{sil} from it and leveraging
#' the special way that the \code{\link[graphics]{plot}} function handles a
#' \code{\link[cluster]{silhouette}} object.
#' 
#' @details This function is intended to be used as part of deciding the ideal
#'   number of clusters to use for analysis. The silhouette score tells
#'   one the goodness of a given clustering looking at a given data points 
#'   dissimilarity within and outside of its cluster. See
#'   \code{\link[cluster]{silhouette}} for further details on its calculation.  
#'   This plot allows one to see how well the individual clusters of a given
#'   clustering perform according to silhouette width. One may pass in the
#'   \code{k_best} element from the \code{multiClust} object for
#'   \code{num_clust}, or use a different value. This is where other
#'   visualizations can be useful.  
#'   See the \emph{Details} section for the \code{\link{multi_clust}} function
#'   and view the \emph{TIP} for a suggested workflow.
#'
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{multiClust} object.
#' @param ... Further arguments to be passed to the \code{\link{plot}} 
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link{multi_clust}}, \code{\link{wss_plot}}, 
#'   \code{\link{gap_plot}}, \code{\link{pca_plot}},
#'   \code{\link{avg_sil_plot}}, \code{\link{clust_boxplot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with sil_plot
#' sil_plot(iris_cluster, num_clust = 3)
sil_plot <- function(clust_obj, num_clust, ...) {
  validate_multi_clust(clust_obj)
  validate_pos_num(list(num_clust = num_clust))
  sil <- clust_obj[["sil"]][[num_clust]]
  plot(sil, main = "Silhouette Plot of K Clusters", ...)
}

#' Plot average silhouette widths for different numbers of clusters.
#' 
#' This function plots the average silhouette widths for all numbers of
#' clusters present in \code{clust_obj}. It provides a wrapper on 
#' \code{\link[graphics]{plot}} by extracting \code{sil_avg} from 
#' \code{clust_obj} and plotting it.
#' 
#' @details This function is intended to be used as part of deciding the ideal
#'   number of clusters to use for analysis. The average silhouette score tells
#'   one the goodness of a given clustering looking at a given data points 
#'   dissimilarity within and outside of its cluster. Typically, one looks for a
#'   peak before which the average silhouette score rises significantly. See
#'   \code{\link[cluster]{silhouette}} for further details on its calculation.
#'   Sometimes, the peak in the score still does not contain enough clusters to
#'   understand the data well. This is where other visualizations can be useful.  
#'   Choosing \code{optimal = TRUE} will circle the optimal number of clusters
#'   based on average silhouette width. See the \emph{Details} section for the
#'   \code{\link{multi_clust}} function and view the \emph{TIP} for a suggested
#'   workflow.
#'
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{clust_gap}.
#' @param optimal Logical. If \code{TRUE}, the optimal number of clusters as
#'   extracted from \code{clust_obj}, based on average silhouette width, is
#'   circled in the plot.
#' @param ... Further arguments to be passed to the \code{\link[graphics]{plot}}
#'   function (besides \code{xlab}, \code{ylab}, \code{main}).
#'
#' @export
#'
#' @seealso \code{\link{multi_clust}}, \code{\link{wss_plot}}, 
#'   \code{\link{gap_plot}}, \code{\link{pca_plot}}, \code{\link{sil_plot}},
#'   \code{\link{clust_boxplot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with avg_sil_plot
#' avg_sil_plot(iris_cluster, optimal = TRUE)
avg_sil_plot <- function(clust_obj, optimal = FALSE, ...) {
  validate_not_null(list(clust_obj = clust_obj, optimal = optimal))
  validate_true_false(list(optimal = optimal))
  validate_multi_clust(clust_obj)
  sil <- clust_obj[["sil_avg"]]
  start <- length(sil[is.na(sil)]) + 1 # We start the krange right after NA's
  end <- length(sil)
  krange <- start:end
  sil <- sil[!is.na(sil)]
  plot(krange, sil, type = "b", xlab = "# of K Clusters",
       ylab = "Average Silhouette Width", 
       main = "Average Silhouette Width by K Clusters",
       ...)
  if (optimal) {
    opti_clust <- clust_obj[["k_best"]]
    index <- opti_clust - (start - 1) # Account for removal of NA's
    points(opti_clust, sil[index], col = "red", pch = 1, cex = 3)
    legend("bottomright", "Optimal Clusters", col = "red", pch = 1)
  }
}

#' Calculate axis label font size for boxplot
#' 
#' This is an internal function that returns a font size for the
#' \code{\link{clust_boxplot}} function that should minimize overlapping of
#' axis labels.
#' @param num_clust The num_clust argument from \code{\link{clust_boxplot}}
#' @return An integer value for the font size
#' @keywords internal
axis_label_size <- function(num_clust) {
  if (num_clust <= 10) {
    return(12)
  } else if (num_clust <= 15) {
    return (10)
  } else if (num_clust <= 20) {
    return (8)
  } 
}

#' Calculate axis label font size for boxplot
#' 
#' This is an internal function that returns a number of columns for the
#' \code{\link{clust_boxplot}} function that should minimize tightness of
#' plot for high number of clusters.
#' @param num_clust The num_clust argument from \code{\link{clust_boxplot}}
#' @return An integer value for the number of columns
#' @keywords internal
boxplot_num_cols <- function(num_clust) {
  if (num_clust < 10) {
    return(5)
  } else if (num_clust <= 15) {
    return (4)
  } else if (num_clust <= 20) {
    return (3)
  } else if (num_clust <= 25) {
    return (2)
  } else if (num_clust <= 30) {
    return (1)
  }
}


#' Plot cluster membership for each feature.
#' 
#' This function plots cluster membership for each feature of \code{d} using box
#' plots. It utlizes \code{facet_wrap} in the \code{\link[ggplot2]{qplot}}
#' function and the desired \code{clust_model} from the argument 
#' \code{clust_obj} to show the box plots of each feature after using the
#' \code{\link[reshape2]{melt}} function to get the correct form of \code{d}.
#' 
#' @details This function is intended to be used to view the shape of a given
#'   clustering to use for analysis. One can view which features (likely
#'   phenotypic markers) are found in each cluster and how the values (likely
#'   relative expression level) are distributed. This function is primarily
#'   intended to try to uncover what sort of underlying structure a given
#'   clustering has found in the data.  
#'   One may pass in the \code{k_best} element from the \code{multiClust} object
#'   for \code{num_clust}, or use a different value. This is where other
#'   visualizations can be useful.  
#'   See the \emph{Details} section for the \code{\link{multi_clust}} function
#'   and view the \emph{TIP} for a suggested workflow.
#'   
#'   \emph{NOTE:} If boolean data is passed in for \code{d}, \code{TRUE} should
#'   be represented with a \code{1} and \code{FALSE} should be represented with
#'   a \code{0}. Also, all columns should be boolean in this case so that the
#'   scaling is correct.
#' 
#' @param d A numeric matrix of data, or an object that can be coerced to
#'   such a matrix (such as a numeric vector or a data frame with all numeric
#'   columns). Note: This should be the same one used to generate 
#'   \code{clust_obj}.
#' @param clust_obj A \code{multiClust} object from which to extract
#'   \code{clust_model} based on the argument \code{num_clust}
#' @param num_clust An integer. The desired number of clusters to be used. Note:
#'   This integer should fall within the krange used to generate the 
#'   \code{multiClust} object.
#' @param ... Further arguments to be passed to the \code{\link{qplot}} 
#'   function (besides \code{xlab}, \code{ylab}).
#'
#' @export
#'
#' @seealso \code{\link{multi_clust}}, \code{\link{wss_plot}}, 
#'   \code{\link{gap_plot}}, \code{\link{pca_plot}}, \code{\link{sil_plot}},
#'   \code{\link{avg_sil_plot}}
#'
#' @examples
#' # First, create a multiClust object
#' library(datasets)
#' iris_cluster <- multi_clust(iris[, 1:4])
#' # Second, use object with clust_boxplot
#' clust_boxplot(iris[, 1:4], iris_cluster, num_clust = 3)
clust_boxplot <- function(d, clust_obj, num_clust, ...) {
  validate_num_data(d)
  validate_multi_clust(clust_obj)
  validate_pos_num(list(num_clust = num_clust))
  meas_vars <- colnames(d)
  d["cluster"] <- clust_obj[["clust_model"]][[num_clust]][["cluster"]]
  if (is_boolean(d[meas_vars])) {
    d_bool <- aggregate(. ~ cluster, data = d, sum)
    d_bool[meas_vars] <- d_bool[meas_vars] / nrow(d)
    m <- reshape2::melt(d_bool, id.vars = "cluster", measure.vars = meas_vars)
    ggplot2::qplot(x = as.factor(cluster), y = value, data = m, geom = "bar", 
                   stat = "identity", fill = as.factor(cluster), xlab = NULL, 
                   ylab = "Proportion expressed", ...) + 
      ggplot2::facet_wrap(~variable, ncol=boxplot_num_cols(num_clust)) + 
      ggplot2::scale_fill_discrete(name = "Cluster") +
      ggplot2::theme(
        axis.text=ggplot2::element_text(size=axis_label_size(num_clust))
      )
  } else {
    m <- reshape2::melt(d, id.vars = "cluster", measure.vars = meas_vars)
    ggplot2::qplot(x = as.factor(cluster), y = value, data = m,
                   geom = "boxplot", fill = as.factor(cluster), xlab = NULL, 
                   ylab = "Relative expression level", ...) + 
      ggplot2::facet_wrap(~variable, ncol=boxplot_num_cols(num_clust)) + 
      ggplot2::scale_fill_discrete(name = "Cluster") +
      ggplot2::theme(
        axis.text=ggplot2::element_text(size=axis_label_size(num_clust))
      )
  }
}
