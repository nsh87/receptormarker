#' @title Get installed version of an R package
#' @description An internal function that retrieves the version of a R package
#' @param pkg The name of the R package
#' @return If the package is installed, returns the version number as a
#'   \code{numeric} (for example, \code{3.8}), otherwise returns NULL.
#' @examples
#' \dontrun{
#' pkg_version('base')
#' }
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
#'   \code{muscle}.
#' @details This is an internal function for checking to see whether or not
#'   the latest version of \code{muscle} is installed. Previous versions
#'   of the package could be installed via CRAN, but the package has now been
#'   removed from CRAN and is only available on Bioconductor. The version
#'   on Bioconductor has different function parameters and therefore users
#'   need to update to the current version in order to prevent the calls
#'   to \code{muscle} from erroring.
#' @return If \code{muscle} is installed the version number will be
#'   returned as a \code{numeric} (for example, \code{3.1}). If the package is
#'   not installed, returns NULL.
#' @keywords internal
muscle_version <- function() {
  pkg_version("muscle")
}


#' @title Check for the existence of Python and Biopython
#' @description Reads the versions of Python and Biopython from this package's
#'   options and warns the user or stops execution depending on the paramater
#'   \code{level}.
#' @param level Either \code{"startup_warn"}, \code{"warn"}, or \code{"stop"},
#'   depending on how you want to handle the non-existence of Python or
#'   Biopython.
#' @return Returns TRUE if the libraries are found, otherwise returns
#'   \code{NULL}.
#' @keywords internal
check_bio_python <- function(level) {
  if (!(level %in% c("startup_warn", "warn", "stop"))) {
    stop("Invalid argument", call.=TRUE)
  }
  
  # Check version numbers store in options
  missing_py_apps <- c()
  if (is.null(getOption("receptormarker.py_version"))) {
    missing_py_apps[length(missing_py_apps) + 1] <- "Python"
  }
  if (is.null(getOption("receptormarker.biopy_version"))) {
    missing_py_apps[length(missing_py_apps) + 1] <- "Biopython"
  }
  
  # Build error message, either singular or plural if both are missing
  if (length(missing_py_apps) > 0) {
    if (length(missing_py_apps) == 1) {
      apps_list <- missing_py_apps
    } else {
      apps_list <- paste0(missing_py_apps, collapse=" and ")
    }
    
    install_biopy <- "http://biopython.org/DIST/docs/install/Installation.html"
    err <- paste0(c("Unable to find",
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
  }
  
  if (exists("err")) {
    if (level == "startup_warn") {
      packageStartupMessage(paste0(c("Warning:", err), collapse=" "))
      return()
    } else if (level == "warn") {
      warning(err, call.=FALSE)
      return()
    } else if (level == "stop") {
      stop(err, call.=FALSE)
    }
  }
  return(TRUE)
}
 

#' @title Check for correct version of MUSCLE
#' @description Reads the version of \code{muscle} from this package's
#'   options and warns the user or stops execution depending on the parameter
#'   \code{level}.
#' @param level Either \code{"startup_warn"}, \code{"warn"}, or \code{"stop"},
#'   depending on how you want to handle the non-existence of the \code{muscle}
#'   package or an outdated version of the package.
#' @return Returns TRUE if the correct version of MUSCLE is found, otherwise
#'   returns \code{NULL}.
#' @keywords internal
check_muscle <- function(level) {
  if (!(level %in% c("startup_warn", "warn", "stop"))) {
    stop("Invalid argument", call.=TRUE)
  }
  
  install_muscle <- paste0(c("http://www.bioconductor.org/packages/release/",
                             "bioc/html/muscle.html"),
                           collapse="")
  if (is.null(getOption("receptormarker.muscle_version"))) {
    err <- paste0(c("The package 'muscle' is not",
                    "installed. Please install the latest",
                    "version from Bioconductor then reload this package. See",
                    install_muscle,
                    "for details."),
                  collapse=" ")
  } else if (getOption("receptormarker.muscle_version") < 3.10 ) {
    err <- paste0(c("The package 'muscle' is not",
                    "at a current version. Please install",
                    "the latest version from Bioconductor then reload this",
                    "package. See", install_muscle, "for details."),
                  collapse=" ")
  }
  
  # Raise message if MUSCLE not found or version is invalid
  if (exists("err")) {
    if (level == "startup_warn") {
      packageStartupMessage(paste0(c("Warning:", err), collapse=" "))
      return()
    } else if (level == "warn") {
      warning(err, call.=FALSE)
      return()
    } else if (level == "stop") {
      stop(err, call.=FALSE)
    }
  }
  return(TRUE)
}


#' @title Get the OS type
#' @description Polls the current system to get the platform type: 'win',
#' 'osx', or 'linux'. The discovery of 'osx' might not be 100%, so if the
#' platform is not discovered to be either one of "win", "osx", but it is found
#' to be "unix" then the return OS type will be "linux". If the platform is
#' not found to be any of the three ("win", "osx", or "unix") then an error
#' is thrown. Sinced OS X is a Unix variation, the default 'linux'
#' should generally be applicable even if the platform is OS X but that
#' couldn't be discovered properly.
#' @return Either 'win', 'osx', or 'linux', as a string.
#' @keywords internal
get_os <- function() {
  os <- "linux"
  if (.Platform$OS.type == "windows") {
    os <- "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    os <- "osx"
  } else if (.Platform$OS.type != "unix") {
    stop("Unknown OS: failed to append HMMER/BLAST+ binaries.", call.=FALSE)
  }
  os
}


#' @title Parse system PATH into a list
#' @description The system \code{PATH} is usually of format
#' \code{/usr/bin:/usr/sbin:/usr/local/bin...} where each folder is separated
#' by a \code{:} character. This function takes the path and splits it into
#' a list of individual paths so that it can be used elsewhere. The function
#' has been taken from the \code{rappdirs} package at:
#' https://github.com/hadley/rappdirs.
#' @return A vector of characters containing the various paths.
#' @keywords internal
parse_path_string <- function(path, sep=":") {
  normalizePath(unique(strsplit(path, sep)[[1]]), mustWork=FALSE)
}


#' @title Get the architecture of the current machine
#' @description Checks to see if a machine uses a x86_64 or ia32 processor.
#' @return Either 'x86_64' or 'ia32', depending on the architecture.
#' @keywords internal
get_architecture <- function() {
  if (sum(grepl("86_64", Sys.info()["machine"])) > 0) {
    machine <-"x86_64"
  } else if (sum(grepl("32", Sys.info()["machine"])) > 0) {
    machine <- "ia32"
  } else {
    err <- paste0(c("Unknown architecture ('x86_64' or 'ia32' sought):",
                    "failed to append HMMER/BLAST+ binaries"),
                  collapse=" ")
    stop(err, call.=FALSE) 
  }
  machine
}


#' @title Add HMMER and NCBI BLAST+ binaries to the system path
#' @description HMMER and BLAST+ binaries are included for multiple platforms.
#' For HMMER, there are 'windows' and 'unix' binaries. For BLAST+, there are
#' 'windows', 'mac-osx', 'linux x86', and 'linux 'linux ia32' binaries. This
#' function detects the current operating system and adds the appropriate
#' binaries to the system \code{PATH} variable if they are not already in
#' the \code{PATH}.
#' @keywords internal
add_binaries_to_path <- function() {
  path_items <- NULL # Will hold vector of items in original system path
  path_sep <- NULL  # Will hold ';' or ':', depending on the operating system
  machine <- NULL  # Will hold 'x86_64' or 'ia32', depending on the architecture
  # Get OS and a list of items in the system path
  os <- get_os()
  if (os == "win") {
    path_items <- parse_path_string(Sys.getenv("PATH"), sep=";")
    path_sep <- ";"
  } else if (os == "osx") {
    path_items <- parse_path_string(Sys.getenv("PATH"))
    path_sep <- ":"
  } else if (os == "linux") {
    path_items <- parse_path_string(Sys.getenv("PATH"))
    path_sep <- ":"
    machine <- get_architecture()  # Holds "x86" or "ia32", needed for BLAST+ 
  }
  
  # Get location of included HMMER binaries in a platform-independent way
  if (os == "win") {
    hmmer <- system.file(file.path("binaries/hmmer/hmmer-3_0-windows/bin"),
                         package="receptormarker")
  } else if (os == "osx") {
    hmmer <- system.file(
      file.path("binaries/hmmer/hmmer-3_1b2-macosx-intel/bin"),
      package="receptormarker")
  } else if (os == "linux" && machine == "x86_64") {
    hmmer <- system.file(
      file.path("binaries/hmmer/hmmer-3_1b2-linux-intel-x86_64/bin"),
      package="receptormarker"
    )
  } else if (os == "linux" && machine == "ia32") {
    hmmer <- system.file(
      file.path("binaries/hmmer/hmmer-3_1b2-linux-intel-ia32/bin"),
      package="receptormarker")
  }
  
  # Get location of included BLAST binaries in a platform independent way
  if (os == "win") {
    blast <- system.file(file.path("binaries/ncbi-blast/windows/bin"),
                         package="receptormarker")
  } else if (os == "osx" || os == "linux") {
    blast <- system.file(file.path("binaries/ncbi-blast/unix/bin"),
                         package="receptormarker")
  } 
  
  # Check if BLAST+ is in path (in a platform independent way) and add it if not
  g <- grepl(file.path("/binaries/ncbi-blast/"), path_items)  # nolint
  if (sum(g) == 0) {
  # Add the BLAST+ binaries to the system path
    Sys.setenv(PATH=paste0(c(Sys.getenv("PATH"), blast), collapse=path_sep)) 
  }
  
  # Check if HMMER is in path (in a platform independent way) and add it if not
  g <- grepl(file.path("/binaries/hmmer-"), path_items)  # nolint
  if (sum(g) == 0) {
  # Add the HMMER binaries to the system path
    Sys.setenv(PATH=paste0(c(Sys.getenv("PATH"), hmmer), collapse=path_sep)) 
  }
}
