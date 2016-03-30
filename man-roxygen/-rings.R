#' @param rings A named character vector that can be used to create outer-ring
#'   annotations on the radial phylogram. There are two ways \code{rings} can
#'   be used.
#'   
#'   1. The names of the vector correspond to column names in the data.frame
#'   \code{d}, and the values should correspond to a desired value in each
#'   column which should be annotated on the ring.
#'   For example, \code{rings=c(FOXP3=1, species="human")} will create
#'   two outer rings, the first of which will be colored whenever the column
#'   "FOXP3" is 1 and the second of which will be colored whenever the column
#'   "species" is "human". Annotations occur on a per-sequence basis when the
#'   PhyloXML represents a non-condensed phylogram. If the PhyloXML represents
#'   a condensed phylogram, annotations occur using individual sequence
#'   populations: if 50\% or more of the cells with a given sequence meet the
#'   current criteria given by \code{rings} then that sequence's ring on the
#'   radial phylogram will be annotated.
#'   
#'   2. Alternatively, \code{rings} can be used to label all unique values in
#'   a single column. In this case, only a single column can be used to create
#'   rings. The name of the vector should correspond to the column name in the
#'   data.frame \code{d}, and the value should be \code{"all"}. For example,
#'   Consider the scenario where column "patient" contains four unique values,
#'   "P1", "P2", "P3", and "P4". Using  \code{rings=c(patient="all")}
#'   will create four outer rings, the first of which will be colored
#'   whenever the column "patient" is "P1", the second whenever it is,
#'   "P2", the third whenever it is "P3", and the fourth whenever it is
#'   "P4". The coloring of the rings is not mutually exclusive, so a
#'   single sequence could have colored rings for, say, "P1" and "P3" at
#'   the same time if there are identical sequences for the first and third
#'   patients in \code{d}. Note that due to the large number of rings this can
#'   create, only a single column can be labeled using this method.
