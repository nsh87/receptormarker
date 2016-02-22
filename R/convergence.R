#' @title Parse the results of the convergence tool into a data frame
#' @description An internal function that extracts the relevant tabular data
#' from the convergence output file and saves it as a data frame for easy
#' manipulation and reading.
#' @param results_file A path to the \code{...-convergence-groups.txt} output
#' file created by the convergence tool. This path is returned by
#' \code{\link{run_convergence}}.
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE} a CSV
#' representation of the returned data frame is written to the
#' \code{verbose_dir}.
#' @template -verbose_dir
#' @return A data frame containing the clusters data. Each row represents a
#' single cluster. The first column contains the number of items in each cluster
#' and each cell in the remaining columns is a sequence belonging to that
#' cluster.
#' @keywords internal
parse_convergence_results <- function(results_file, verbose, verbose_dir) {
  # Since there are no headers, you can't parse the CSV and have blank cells
  # where you need them. So you need to figure out how many columns the data
  # frame should have, then create a data frame and specify that many columns.
  d <- read.csv(results_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  largest_cluster_size <- max(d[, 1])
  num_cols <- largest_cluster_size + 2  # Need to account for 1st and 2nd cols
  d <- read.csv(results_file, sep="", header=FALSE,
                col.names=as.character(c(1:num_cols)), stringsAsFactors=FALSE)
  
  # Write data frame to verbose dir if user wants it
  if (verbose) {
    csv_path <- gsub(".txt$", "-parsed.csv", results_file)
    write.table(d, csv_path, sep=",", row.names=FALSE)
    # CSV gets written to convergence's temporary dir, need to copy it
    if (file.exists(csv_path)) {
      file.copy(csv_path, verbose_dir)
    }
  }
  d
}


#' @title Create an interactive clustering of adaptive repertoire convergence
#' @description An internal function that utilizes Jacob Glanville's
#' \code{convergence.pl} script and associated files to find statistical
#' clustering of adaptive repertoire convergence.
#' @details Requires HMMER and blastn libraries. These binaries are included
#' in the R package, but may present difficulties on non-Unix computers.
#' @template -seqs
#' @param verbose \code{TRUE} or \code{FALSE}, depending on whether or not the
#' output of the convergence script should be printed and the output file should
#' be copied to the \code{verbose_dir}.
#' @template -verbose_dir
#' @return A path to the output file containing the convergence groups in
#' tabular format.
#' @keywords internal
run_convergence <- function(seqs, verbose, verbose_dir){
  seqs <- unique(seqs)
  
  # Write sequences to file and create temporary dir to hold input/output files
  convergence_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(convergence_tmpdir)
  seqs_file <- tempfile(pattern="sequences-deduped-", tmpdir=convergence_tmpdir,
                        fileext=".txt")
  write(seqs, seqs_file)
  # Copy deduped sequences file to output dir if user wants it
  if (verbose && file.exists(seqs_file)) {
    file.copy(seqs_file, verbose_dir)
  }
  
  # Add the HMMER and blastn binaries to the system path
  original_sys_path <- Sys.getenv("PATH")
  hmmer <- system.file("binaries/hmmer-3.1b2/binaries",
                       package = "receptormarker")
  blastn <- system.file("binaries/ncbi/blast/bin", package = "receptormarker")
  Sys.setenv(PATH=paste0(c(original_sys_path, hmmer, blastn), collapse=":"))
  
  # Run the convergence tool
  tryCatch({
    convergence <- system.file("perl/vdjfasta/bin", "convergence-pipeline.pl",
                               package="receptormarker")
    system(sprintf("perl %s --textfile=%s", convergence, seqs_file),
           ignore.stdout=!verbose, ignore.stderr=!verbose)
    # Collect output
    output_file <- gsub(".txt$", "-convergence-groups.txt", seqs_file)
    # Copy output file to verbose dir if user wants it
    if (verbose && file.exists(seqs_file)) {
      file.copy(output_file, verbose_dir)
    }
    # Revert system path back to what it was originally 
    Sys.setenv(PATH=original_sys_path)
  },
  error = function(e) {
    # Need to set system path back to what it was if something goes wrong
    Sys.setenv(PATH=original_sys_path)
    stop(e, call.=FALSE) 
  }
  )
  output_file
}

#' Group adaptive repertoire convergence by paratope hotspots for CDRs
#' @description Generates convergence motifs and establishes a statistical
#' cutoff to determine convergence groups. Convergence groups can be plotted
#' using \code{\link{convergence_plot}}.
#' @details Inspection of the convergence groups can (and likely should) be
#' performed by analyzing the data frame contained within the returned object.
#' See the examples below for more info.
#' @template -d
#' @template -seqs_col
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE}, additional output
#' is printed to the console and the sequences and resulting cluster
#' (convergence groups) are written to a folder in the working directory.
#' @return An object of class \code{\link{convergenceGroups-class}} that can be 
#' inspected or used to generate network graphs using
#' \code{\link{convergence_plot}}.
#' @examples 
#' data(tcr)  # Packaged data set, a data.frame from a CSV file
#' tcr_reduced <- tcr[1:100, ]
#' converged <- convergence(tcr_reduced, seqs_col='seqs')
#' # You can then: View(converged@groups)
#' @export
convergence <- function(d, seqs_col=NULL, verbose=FALSE) {
  validate_not_null(list(d=d, seqs_col=seqs_col, verbose=verbose))
  validate_d_seqs(d, seqs_col)
  
  # Create verbose dir
  if (verbose) {
    verbose_dir <- tempfile("convergence-", tmpdir=getwd(), fileext="")
    dir.create(verbose_dir)
  }
  
  # Step 1: Clean the user-supplied data and the sequences
  clean <- clean_data(d, seqs_col, verbose, verbose_dir, verbose_format="txt")
  seqs <- clean[["seqs"]]
  
  # Step 2: Run the convergence tool
  results_file <- run_convergence(seqs, verbose, verbose_dir)
  
  # Step 3: Parse output file into a data frame
  groups <- parse_convergence_results(results_file, verbose, verbose_dir)
  
  # Step 4: Rename columns c("num_items", "group_name", "X1", "X2", "X3", etc.)
  num_cols <- ncol(groups)
  if (num_cols >= 2) {
    colnames(groups)[1:2] <- c("num_items", "group_name")
  }
  if (num_cols >= 3) {
    num_group_cols <- length(3:num_cols)
    numbered_cols <- lapply(1:num_group_cols, function(x) {
      paste0("X", x, sep="") 
    })
    colnames(groups)[3:num_cols] <- numbered_cols
  }
  
  # Step 5: Sort the convergence groups from largest to smallest
  num_rows <- nrow(groups)
  if (num_rows > 0) {
    groups <- groups[order(-groups["num_items"]), ]
    row.names(groups) <- 1:num_rows
  }
    
  # Instantiate a 'convergenceClust' class to hold clusters and return it
  new("convergenceGroups", groups=groups)
}
