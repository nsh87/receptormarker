#' Analyze antibody receptor sequences and phenotypic markers.
#'
#' Useful for quick normalization, clustering, and creating phylogenetic trees.
#' @docType package
#' @name receptormarker
NULL


#' @title Creates and stores package options using R's \code{.onLoad} hook.
#' @description Saved options will be available throughout the package using
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
    receptormarker.py_version = py_version(),
    receptormarker.biopy_version = biopy_version(), 
    receptormarker.muscle_version = muscle_version()
  )
  toset <- !(names(op.receptormarker) %in% names(op))
  if(any(toset)) options(op.receptormarker[toset])
  
  invisible()
}


# Check if Python/R libraries are installed.
# For example, presents user with warning if Biopython, Python, or muscle is not
# installed, but allows user to continue using the package.
.onAttach <- function(libname, pkgname) {
  check_bio_python(level="startup_warn")
  check_muscle(level="startup_warn")
}
