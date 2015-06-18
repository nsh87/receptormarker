#' @title T-cell data set
#' @description A data set containing T-cell receptor sequences and
#' phenotypic markers. Each row corresponds to data for a single cell.
#' @format A data frame with 265 rows and 23 columns:
#' \describe{
#'   \item{stim}{whether or not the cell received stimulation}
#'   \item{seqs}{T-cell receptor (TCR) sequences}
#'   \item{FOXP3, GATA3, etc.}{phenotypic markers - the absolute measurements
#'     were thresholded in order to make all markers either NA or 1. Presence
#'     of the marker is indicated by 1.}
#' }
#' @usage data(tcr)
#' @references Han, A., Glanville, J., Hansmann, L., & Davis, M. M. (2014).
#' Linking T-cell receptor sequence to functional phenotype at the single-cell
#' level. \emph{Nature biotechnology}.
#' @keywords datasets
#' @docType data
"tcr"
