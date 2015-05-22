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
    receptormarker.py_version = py_version(),  
    receptormarker.biopy_version = biopy_version()
  )
  toset <- !(names(op.receptormarker) %in% names(op))
  if(any(toset)) options(op.receptormarker[toset])
  
  invisible()
}


#' \code{library(receptormarker)} and presents user with warning if none found.
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  missing_apps = c()
  if (is.null(getOption('receptormarker.py_version'))) {
    missing_apps[length(missing_apps) + 1] <- "Python"
  }
  if (is.null(getOption('receptormarker.biopy_version'))) {
    missing_apps[length(missing_apps) + 1] <- "Biopython"
  }
  
  if (length(missing_apps) > 0) {
    if (length(missing_apps) == 1) {
      apps_list <- missing_apps
    } else {
      apps_list <- paste0(missing_apps, collapse=" and ")
    }
    
    install_biopy <- 'http://biopython.org/DIST/docs/install/Installation.html'
    missing_apps_warning <- paste0(c("Warning: Unable to find",
                                     apps_list,
                                     "on your system. In order to achieve the",
                                     "best results it is strongly suggest to",
                                     "install Python's Biopython to your",
                                     "system's path. See",
                                     install_biopy,
                                     "for installation instructions. Once",
                                     "installed, please reload this package."
                                     ),
                                   collapse=" "
    )
    packageStartupMessage(missing_apps_warning)
  }
}