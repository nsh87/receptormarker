#' @title Create a path to a tempory XML file for plotting with Cytoscape
#' @description An internal function that creates a temporary directory to hold
#' the XML file for a Cytoscape network diagram. The directory created is only
#' intended to hold the XML, and no other files, because in order for
#' htmlwidgets to make the XML available in the HTML it copies the entire parent
#' directory.
#' @return A filepath that can be used to save a XML.
#' @keywords internal
convergence_xml_path <- function() {
  # Create a temp dir for the XML file
  convergence_xml_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(convergence_xml_tmpdir)  # htmlwidgets copies entire dir to browser
  
  xml_file <- tempfile(pattern="phyloxml-", tmpdir=convergence_xml_tmpdir,
                       fileext=".xml")
}


#' @title Create a XML file to plot with Cytoscape
#' @description An internal function that parses the output of the convergence
#' tool and creates an XML file to represent a specific row (cluster) in the
#' file.
#' @return A path to the Cytoscape XML file.
#' @keywords internal
cytoscape_xml <- function(d, row_num, verbose, verbose_dir) {
  xml_file <- convergence_xml_path()
  
  # Get cluster nodes, label, and size
  num_cols <- ncol(d)
  cluster <- d[row_num, c(3:num_cols)]
  nodes <- na.omit(cluster[cluster != ''])
  cluster_label <- d[row_num, 2]
  num_nodes <- d[row_num, 1]
  # Ensure you have the correct number of nodes
  if (num_nodes != length(nodes)) {
    stop("Number of nodes parsed does not match expected number of nodes.",
         call.=FALSE)
  }
  
  # Create the XML tree
  root <- XML::newXMLNode("graphml")
  label_key <- XML::newXMLNode("key", 
                               attrs=c(id="label", "for"="all",
                                       "attr.name"="label",
                                       "atr.type"="string"),
                               parent=root)
  weight_key <- XML::newXMLNode("key",
                                attrs=c(id="weight", "for"="node",
                                        "attr.name"="weight",
                                        "attr.type"="double"),
                                parent=root)
  graph_node <- XML::newXMLNode("graph",
                                attrs=c(id="0", "edgedefault"="undirected",
                                        "label"=cluster_label),
                                parent=root)
  # Add cluster nodes into the XML tree as children of graph_node, <graph>
  lapply(c(1:num_nodes), function(x) {
    single_node <- XML::newXMLNode("node", attrs=c(id=x), parent=graph_node)
    XML::newXMLNode("data", nodes[[x]], attrs=c("key"="label"),
                    parent=single_node)
  })
  
  # Add edges to XML, connecting each node to every other node
  lapply(c(1:(num_nodes - 1)), function(x) {
    lapply(c((x + 1):num_nodes), function(y) {
      XML::newXMLNode("edge", attrs=c("source"=x, "target"=y),
                      parent=graph_node)
    })
  })
  
  # Write XML to file
  XML::saveXML(root, file=xml_file,
          prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
          indent=FALSE)
  
  # Copy XML file to verbose dir if user wants it
  if (verbose && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  
  xml_file
}


#' @title Parse the results of the convergence tool into a data frame
#' @description An internal function that extracts the relevant tabular data
#' from the convergence output file and saves it as a data frame for easy
#' manipulation and reading.
#' @param results_file A path to the \code{...-convergence-groups.txt} output
#' file created by the convergence tool. This path is returned by
#' \code{\link{run_convergence}}.
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


#' Create an interactive convergence network diagram
#'
#' Creates a JavaScript-based network diagram in RStudio or a browser window.
#' This visualization can help identify clusters of functionally similar
#' complementary determining regions (CDRs) of antibody or T-cell receptors.
#'  
#'
#' @import htmlwidgets
#' @example 
#' 
#' @export
convergence <- function(d, seqs_col=NULL, browser=FALSE, verbose=FALSE) {
  validate_not_null(list(d=d, browser=browser, verbose=verbose))  
  validate_d_seqs(d, seqs_col)
  
  # Not necessary to have as func parameters; these will get set automatically
  width <- NULL
  height <- NULL
  
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
  clusters <- parse_convergence_results(results_file, verbose, verbose_dir)
  
  # Step 4: Save a specific cluster (row) to XML to send to Cytoscape
  # Plot the largest cluster, or the first row if all clusters are the same size
  largest_cluster <- which.max(clusters[, 1])
  xml_file <- cytoscape_xml(clusters, row_num=largest_cluster, verbose,
                            verbose_dir)
  
  # Add the XML file as an HTML dependency so it can get loaded in the browser
  convergence_xml <- htmltools::htmlDependency(
    name = "convergence_xml",
    version = "1.0",
    src = c(file=dirname(xml_file)),
    attachment = list(xml=basename(xml_file))
  )
  
  # Create widget
  htmlwidgets::createWidget(
    name = "convergence",
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(
      padding = 22,
      viewer.suppress = browser,
      browser.fill = TRUE
    ),
    package = "receptormarker",
    dependencies = convergence_xml
  )
}

#' Widget output function for use in Shiny
#'
#' @export
convergenceOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'convergence', width, height,
                    package = 'receptormarker')
}

#' Widget render function for use in Shiny
#'
#' @export
renderConvergence <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, convergenceOutput, env, quoted = TRUE)
}
