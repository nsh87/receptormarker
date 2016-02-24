#' @title Clean data for a phylogram or convergence and extract the sequences
#' @description An internal function that removes rows of data where the
#' sequence column is \code{NA} or an empty string.
#' @template -d
#' @template -seqs_col
#' @param verbose \code{TRUE} or \code{FALSE}, depending on whether or not the
#' cleaned sequences should be written to the \code{verbose_dir}.
#' @template -verbose_dir
#' @param verbose_format \code{"fasta"} or \code{"txt"}, indicating how the
#' cleaned sequences should be written to the \code{verbose_dir}. The options
#' correspond to FASTA format and a simple text file with each sequence on
#' its own line. If \code{verbose} is \code{FALSE}, nothing will be written.
#' @return A named list containing the cleaned \code{"seqs"} and \code{"d"}.
#' @keywords internal
clean_data <- function(d, seqs_col, verbose, verbose_dir, verbose_format) {
  # If d is a >1D data frame
  if (class(d) == "data.frame" && dim(d)[2] > 1) {
    d_clean <- d[d[, seqs_col] != "", ]  # Remove rows with no sequences
    d_clean <- d_clean[complete.cases(d_clean[, seqs_col]), ]  # Remove NA rows
    d_clean <- d_clean[with(d_clean, order(d_clean[, seqs_col])), ]
    seqs <- as.character(d_clean[, seqs_col])
    # If d is a 1D data frame or vector
  } else if ((class(d) == "data.frame" && dim(d)[2] == 1) ||
             (class(d) == "character")) {
    d_clean <- d[d != ""]  # Remove blank sequences
    d_clean <- d_clean[!is.na(d_clean)]  # Remove NA's
    d_clean <- sort(d_clean)
    seqs <- as.character(d_clean)
  }
  # Write the sequences if the user wants them
  if (verbose) {
    if (class(d_clean) == "data.frame") {
      if (verbose_format == "fasta") {
        seqs_output <- with(d_clean,
                            paste0(">", seqs, "\n", seqs, collapse="\n"))
      } else if (verbose_format == "txt") {
        seqs_output <- with(d_clean, paste0(seqs, collapse="\n"))
      }
    } else if (class(d_clean) == "character") {
      if (verbose_format == "fasta") {
        seqs_output <- paste0(">", seqs, "\n", seqs, collapse="\n")
      } else if (verbose_format == "txt") {
        seqs_output <- paste0(seqs, collapse="\n")
      }
    }
    seqs_file <- tempfile(pattern="sequences-", tmpdir=verbose_dir,
                          fileext=".txt")
    write(seqs_output, seqs_file)
  }
  list("seqs"=seqs, "d"=d_clean)
}
