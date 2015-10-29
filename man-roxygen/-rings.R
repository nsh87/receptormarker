#' @param rings A named character vector that can be used to create outer-ring
#'   annotations on the radial phylogram. The names of the vector must correspond
#'   to column names in the data.frame \code{d}, and the values should
#'   correspond to a desired value in each column which should be annotated
#'   on the ring. For example, \code{c(FOXP3=1, species="human")} will create
#'   two outer rings, the first of which will be colored whenever the column
#'   "FOXP3" is 1 and the second of which will be colored whenever the column
#'   "species" is "human". Annotations occur on a per-sequence basis when the
#'   PhyloXML represents a non-condensed phylogram. If the PhyloXML represents
#'   a condensed phylogram, annotations occur using individual sequence
#'   populations: if 50\% or more of the cells with a given sequence meet the
#'   current criteria given by \code{rings} then that sequence's ring on the
#'   radial phylogram will be annotated.
