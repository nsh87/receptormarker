#' Analyze antibody receptor sequences and phenotypic markers.
#'
#' Useful for quick normalization, clustering, and creating phylogenetic trees.
#'
#' @docType package
#' @name receptormarker
NULL

#' Creates and stores package options using R's \code{.onLoad} hook.
#' 
#' Saved options will be available throughout the package using
#' \code{getOptions(receptormarker.option_name)}. Options include:
#' \enumerate{
#'   \item \code{py_version}: Python's version; empty string if not installed.
#'     Useful for checking if Python is installed during \code{.onAttach()}
#'     hook when package is loaded using \code{library()}.
#'   \item \code{biopy_version}: Biopython's version; emptry string if not
#'     installed. Useful for checking if Biopython is installed in any functions
#'     that use the Biopython package.
#' }
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.receptormarker <- list(
    devtools.py_version = py_version(),  
    devtools.biopy_version = biopy_version()
  )
}