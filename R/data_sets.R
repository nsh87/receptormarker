#' @title T-cell data set
#' @description A data set containing T-cell receptor sequences and
#' phenotypic markers. Each row corresponds to data for a single cell.
#' @details A threshold has been applied to the absolute measurements for the
#' phenotypic markers in order to make their values `NA` or 1. Presence of a
#' marker is indicated by 1.
#' @format A data frame with 265 rows and 23 columns:
#' \describe{
#'   \item{group}{whether or not the cell received stimulation}
#'   \item{seqs}{T-cell receptor (TCR) sequences}
#'   \item{FOXP3}{phenotypic marker}
#'   \item{GATA3}{phenotypic marker}
#'   \item{GZMB}{phenotypic marker}
#'   \item{IFNG}{phenotypic marker}
#'   \item{IL10}{phenotypic marker}
#'   \item{IL12A}{phenotypic marker}
#'   \item{IL13}{phenotypic marker}
#'   \item{IL17A}{phenotypic marker}
#'   \item{IL2}{phenotypic marker}
#'   \item{IL21}{phenotypic marker}
#'   \item{IL4}{phenotypic marker}
#'   \item{IL5}{phenotypic marker}
#'   \item{IL9}{phenotypic marker}
#'   \item{PRF1}{phenotypic marker}
#'   \item{RORC}{phenotypic marker}
#'   \item{TBET}{phenotypic marker}
#'   \item{TGFB1}{phenotypic marker}
#'   \item{TNF}{phenotypic marker}
#'   \item{BCL6}{phenotypic marker}
#'   \item{RUNX1}{phenotypic marker}
#'   \item{RUNX3}{phenotypic marker}
#' }
#' @usage data(tcr)
#' @references Han, A., Glanville, J., Hansmann, L., & Davis, M. M. (2014).
#' Linking T-cell receptor sequence to functional phenotype at the single-cell
#' level. \emph{Nature biotechnology}.
#' @keywords datasets
#' @docType data
"tcr"
