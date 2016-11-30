#' @title Parse the groups output by the convergence tool into a data frame
#' @description An internal function that extracts the relevant tabular data
#' from the convergence output file and saves it as a data frame for easy
#' manipulation and reading. Creates meaningful column names and sorts the
#' groups from largest to smallest.
#' @param groups_file A path to the \code{...-convergence-groups.txt} output
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
parse_convergence_groups <- function(groups_file, verbose, verbose_dir) {
  # Since there are no headers, you can't parse the CSV and have blank cells
  # where you need them. So you need to figure out how many columns the data
  # frame should have, then create a data frame and specify that many columns.
  d <- read.csv(groups_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  largest_cluster_size <- max(d[, 1])
  num_cols <- largest_cluster_size + 2  # Need to account for 1st and 2nd cols
  d <- read.csv(groups_file, sep="", header=FALSE,
                col.names=as.character(c(1:num_cols)), stringsAsFactors=FALSE)
  
  # Rename columns c("num_items", "group_name", "X1", "X2", etc.)
  num_cols <- ncol(d)
  if (num_cols >= 2) {
    colnames(d)[1:2] <- c("num_items", "group_name")
  }
  if (num_cols >= 3) {
    num_group_cols <- length(3:num_cols)
    numbered_cols <- lapply(1:num_group_cols, function(x) {
      paste0("X", x, sep="") 
    })
    colnames(d)[3:num_cols] <- numbered_cols
  }
  
  # Sort the convergence groups from largest to smallest
  num_rows <- nrow(d)
  if (num_rows > 0) {
    d <- d[order(-d["num_items"]), ]
    row.names(d) <- 1:num_rows
  }
  
  # Write data frame to verbose dir if user wants it
  if (verbose) {
    csv_path <- gsub(".txt$", "-parsed.csv", groups_file)
    write.table(d, csv_path, sep=",", row.names=FALSE)
    # CSV gets written to convergence's temporary dir, need to copy it
    if (file.exists(csv_path)) {
      file.copy(csv_path, verbose_dir)
    }
  }
  d
}


#' @title Parse the network output by the convergence tool into a data frame
#' @description An internal function that puts the network file output by the
#' tool into a data frame. The network file contains information about which
#' nodes should connect to which other nodes.
#' @details The first and second columns of the data capture which nodes
#' connect to each other. If there is a node in a group (not this data, but in
#' the other file output by convergence) it should be looked up in this network
#' data to see what else the node connects to. Some nodes have no connections,
#' i.e. if they are not found in this network information. The 3rd column in
#' the network info contains either 'local', 'global', or 'singleton'.
#' Singletons have no connections, and local should just be a different color
#' than global connections in the convergence graph.
#' @param network_file A path to the \code{...-clone-network.txt} output
#' file created by the convergence tool. This path is returned by
#' \code{\link{run_convergence}}.
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE} a CSV
#' representation of the returned data frame is written to the
#' \code{verbose_dir}.
#' @template -verbose_dir
#' @return A data frame containg the network data.
#' @keywords internal
parse_convergence_network <- function(network_file, verbose, verbose_dir) {
  d <- read.csv(network_file, sep="", header=FALSE, stringsAsFactors=FALSE,
                col.names=c("node1", "node2", "type"))
  
  # Write data frame to verbose dir if user wants it
  if (verbose) {
    csv_path <- gsub(".txt$", "-parsed.csv", network_file)
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
#' in the R package, but may present difficulties on non-Unix computers. These
#' binaries are added to the system \code{PATH} when the package is loaded.
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
  write(seqs, seqs_file, sep="\n")
  # Copy deduped sequences file to output dir if user wants it
  if (verbose && file.exists(seqs_file)) {
    file.copy(seqs_file, verbose_dir)
  }
  
  # Run the convergence tool
  convergence <- system.file(file.path("perl", "vdjfasta", "bin"),
                             "gliph-group-discovery.pl",
                             package="receptormarker", mustWork=TRUE)
  cmd <- sprintf("perl %s --textfile=%s", convergence, seqs_file)
  system(cmd, ignore.stdout=!verbose, ignore.stderr=!verbose)
  
  # Collect output file containing groups
  groups_file <- gsub(".txt$", "-convergence-groups.txt", seqs_file)
  if (!file.exists(groups_file)) {
    stop("Convergence output not found: no convergence groups, or error.",
         call.=FALSE)
  }
  
  # Collect output file containing network
  network_file <- gsub(".txt$", "-clone-network.txt", seqs_file)
  if (!file.exists(network_file)) {
    stop("Convergence network not found: no node link, or error.",
         call.=FALSE)
  }
  
  # Copy output file to verbose dir if user wants it
  if (verbose && file.exists(groups_file) && file.exists(network_file)) {
    file.copy(groups_file, verbose_dir)
    file.copy(network_file, verbose_dir)
  }
  
  list("groups_file"=groups_file, "network_file"=network_file)
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
  results <- run_convergence(seqs, verbose, verbose_dir)
  groups_file <- results[["groups_file"]]
  network_file <- results[["network_file"]]
  
  # Step 3: Parse convergence output groups file into a data frame
  groups <- parse_convergence_groups(groups_file, verbose, verbose_dir)
  
  # Step 4: Parse the network information output file into a data frame
  network <- parse_convergence_network(network_file, verbose, verbose_dir)
  
  # Instantiate a 'convergenceClust' class to hold clusters and return it
  new("convergenceGroups", groups=groups, network=network)
}
