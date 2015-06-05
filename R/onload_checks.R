#' @title Get installed version of an R package
#' @description An internal function that retrieves the version of a R package
#' @param pkg The name of the R package
#' @return If the package is installed, returns the version number as a
#'   \code{numeric} (for example, \code{3.8}), otherwise returns NULL.
#' @examples
#'   pkg_version('base')
#' @keywords internal
pkg_version <- function(pkg) {
  tryCatch({
    v <- as.character(packageVersion(pkg))
    vrsn <- strsplit(v, ".", fixed=TRUE)[[1]][1:2]
    vrsn <- paste0(vrsn, collapse=".")
    vrsn <- as.numeric(vrsn)
  },
  error = function(e) {
    return(NULL)
  }
  )
}


#' @title Get installed version of external application
#' @description  An internal function that retrieves the version of an external
#'   application by executing code on the system's command line.
#' @param cmd The system command to run, as a single string. The command
#'   must generate output in the form of \emph{"Application version"}, such
#'   as \code{Python 2.7.6}. The system command must be a single line, and
#'   therefore should use characters such as \code{;} to indicate the end of
#'   lines if necessary in the command.
#' @param application The name of an application, as a string, to check for in
#'   the output of the system command. If, for example, "Python" is
#'   given, the version number will only be returned if the output of the system
#'   command begins with \code{Python}.
#' @details If necessary, the command passed to the paramater \code{cmd} should
#'   manipulate the system's output so that it matches the output format
#'   \emph{"Application version"}.
#'   
#'   Uses the \code{\link[base]{system}} function to execute system commands.
#'   WARNING: The \code{system} function does not always capture output to the
#'   Terminal window into a variable. For example, running
#'   
#'   \code{resp <- system(sprintf('python --version'), intern=TRUE)}
#'   
#'   will not capture the version number to the variable \code{resp}, but
#'   invoking a Python command will (see Examples below). Before using this
#'   function make sure you can capture the output of \code{system} to a
#'   variable using your \code{cmd}.
#' @return If the package is installed, returns the version number as a  
#'   \code{numeric} (for example, \code{2.7}), otherwise returns NULL.
#' @examples
#' \dontrun{
#' cmd <- paste0('python -c "import platform;',
#'               'print(\'Python {}\'.format(platform.python_version()))"',
#'                collapse=' ')
#' versn(cmd, 'Python')
#' }
#' @keywords internal
versn <- function(cmd, application) {
  suppressWarnings(
    sys_cmd <- system(sprintf(cmd), intern=TRUE, ignore.stderr=TRUE)
  )
  if (length(sys_cmd) == 0) {
    return()
  }
  tryCatch({
      response <- strsplit(sys_cmd, " ")[[1]]
      if (response[1] == application) {
        full_version <- strsplit(response[2], "\\.")[[1]]  # Uses regex
        versn <- as.numeric(paste0(full_version[c(1:2)], collapse="."))
      } else {
        return()
      }
    },
    warning = function(w) {
      return()
    },
    error = function(e) {
      return()
    }
  )
  return(versn)
}

#' @title Get Python version
#' @description Get the major and minor version of the system's Python
#'   distribution.
#' @details This is an internal function for checking to see whether or not
#'   Python is installed on the system during the \code{.onLoad()} hook.
#' @return If Python is installed the version number will be returned as a
#'   \code{numeric} (for example, \code{2.7}). If Python is not installed,
#'   returns \code{NULL}.
#' @keywords internal
py_version <- function() {
  cmd <- paste0("python -c \"import platform;",
                "print('Python {}'.format(platform.python_version()))\"",
                collapse=" ")
  versn(cmd, "Python")
}


#' @title Get Biopython version
#' @description Get the major and minor version of the installed Biopython
#'   package.
#' @details This is an internal function for checking to see whether or not
#'   Biopython is installed on the system during the \code{.onLoad()} hook.
#' @return If Biopython is installed the version number will be returned as a
#'   \code{numeric} (for example, \code{1.62}). If Biopython is not installed,
#'   returns \code{NULL}.
#' @keywords internal
biopy_version <- function() {
  command <- paste0("python -c \"import Bio;",
                "print('Biopython {}').format(Bio.__version__)\"",
                collapse=" ")
  versn(command, "Biopython") 
}


#' @title Get MUSCLE version
#' @description Get the major and minor version of the R package
#'   \code{\link{muscle}}.
#' @details This is an internal function for checking to see whether or not
#'   the latest version of \code{\link{muscle}} is installed. Previous versions
#'   of the package could be installed via CRAN, but the package has now been
#'   removed from CRAN and is only available on Bioconductor. The version
#'   on Bioconductor has different function parameters and therefore users
#'   need to update to the current version in order to prevent the calls
#'   to \code{\link{muscle}} from erroring.
#' @return If \code{\link{muscle}} is installed the version number will be
#'   returned as a \code{numeric} (for example, \code{3.1}). If the package is
#'   not installed, returns NULL.
muscle_version <- function() {
  pkg_version('muscle')
}
