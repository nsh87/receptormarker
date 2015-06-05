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


#' @title Check if Biopython is installed on \code{library(receptormarker)}
#' @description Presents user with warning if Biopython or Python is not
#'   installed, but allows user to continue using the package.
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Print warning about missing Python applications
  missing_py_apps <- c()
  if (is.null(getOption("receptormarker.py_version"))) {
    missing_py_apps[length(missing_py_apps) + 1] <- "Python"
  }
  if (is.null(getOption("receptormarker.biopy_version"))) {
    missing_py_apps[length(missing_py_apps) + 1] <- "Biopython"
  }
  
  if (length(missing_py_apps) > 0) {
    if (length(missing_py_apps) == 1) {
      apps_list <- missing_py_apps
    } else {
      apps_list <- paste0(missing_py_apps, collapse=" and ")
    }
    
    install_biopy <- "http://biopython.org/DIST/docs/install/Installation.html"
    missing_py_apps_warning <- paste0(c("Warning: Unable to find",
                                     apps_list,
                                     "on your system. In order to achieve the",
                                     "best results it is strongly suggested to",
                                     "install Python's Biopython to your",
                                     "system's path. See",
                                     install_biopy,
                                     "for installation instructions. Once",
                                     "installed, please reload this package.",
                                     "If you use Biopython through Enthought",
                                     "Canopy or Anaconda, these tools run in",
                                     "an environment on top of your system",
                                     "path, which is not accessible through",
                                     "R; please install Biopython to your",
                                     "system path, as well."
                                     ),
                                   collapse=" ")
    packageStartupMessage(missing_py_apps_warning)
  }
  
  # Print warning about any missing R packages from Bioconductor
  install_muscle <- paste0(c("http://www.bioconductor.org/packages/release/",
                             "bioc/html/muscle.html"),
                           collapse="")
  if (is.null(getOption("receptormarker.muscle_version"))) {
    missing_r_apps_warning <- paste0(c("Warning: The package 'muscle' is not",
                                       "installed. Please install the latest",
                                       "version from Bioconductor. See",
                                       install_muscle,
                                       "for details."),
                                     collapse=" ")
    packageStartupMessage(missing_r_apps_warning)
  } else if (getOption("receptormarker.muscle_version") < 3.10 ) {
    missing_r_apps_warning <- paste0(c("Warning: The package 'muscle' is not",
                                       "at a current version. Please install",
                                       "the latest version from Bioconductor.",
                                       "See", install_muscle, "for details."),
                                     collapse=" ")
    packageStartupMessage(missing_r_apps_warning)
  }
}
