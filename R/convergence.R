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
  
  xml_file <- tempfile(pattern="convergence_xml-",
                       tmpdir=convergence_xml_tmpdir,
                       fileext=".xml")
}


#' @title Create a XML file to plot with Cytoscape
#' @description An internal function that parses the output of the convergence
#' tool and creates an XML file to represent a specific row (cluster) in the
#' file.
#' @param d A data frame of clusters, typically returned by
#' \code{\link{run_convergence}}.
#' @param row_num An integer indicating which row number in \code{d} to graph.
#' Each row should correspond to a single cluster.
#' @param labels \code{TRUE} or \code{FALSE}, depending on whether or not node
#' labels should be written to the XML and therefore displayed in Cytoscape.
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE} the XML file is
#' coped to the \code{verbose_dir}.
#' @template -verbose_dir
#' @return A path to the Cytoscape XML file.
#' @keywords internal
cytoscape_xml <- function(d, row_num, labels, verbose, verbose_dir) {
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
  schema_location <- paste0(
    c("http://graphml.graphdrawing.org/xmlns",
      "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"
      ), collapse=" "
  )
  root <- XML::newXMLNode(
    "graphml",
    namespaceDefinitions=c(
      "http://graphml.graphdrawing.org/xmlns",
      "xsi"="http://www.w3.org/2001/XMLSchema-instance"
    ),
    attrs=c("xsi:schemaLocation"=schema_location)
  )
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
  # If 'labels' is FALSE, don't show labels by using empty string for labels
  if (!labels) {
    nodes <- rep("", length(nodes))
  }
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


#' @title Create an interactive convergence network diagram of similar CDRs
#' @description Creates a JavaScript-based network diagram in RStudio or a
#' browser window. This visualization can help identify clusters of
#' functionally similar complementary determining regions (CDRs) of antibody or
#' T-cell receptors.
#' @template -d
#" @template -seqs_col
#' @param background_color A valid hex color code as a string for the canvas
#' color.
#' @param node_shape A string indicating what shape to use for the nodes. One of
#' \code{"ELLIPSE"}, \code{"RECTANGLE"}, \code{"TRIANGLE"}, \code{"DIAMOND"},
#' \code{"HEXAGON"}, \code{"OCTAGON"}, \code{"PARALLELOGRAM"},
#' \code{"ROUNDRECT"}, \code{"VEE"}.
#' @param border_width An integer border width for the nodes.
#' @param border_color A valid hex color code as a string for the border color
#' of the nodes.
#' @param node_color A valid hex color code as a string for the node color.
#' @param node_size An integer representing node size.
#' @param labels A logical indicating whether or not node labels should be
#' displayed. If \code{TRUE}, labels will be the sequences used for analysis in
#' the convergence pipeline.
#' @param label_vertical_pos A string indicating the vertical position of node
#' labels. Must be one of \code{"top"}, \code{"middle"}, or \code{"bottom"}.
#' Ignored if \code{label=FALSE}.
#' @param label_horizontal_pos  A string indicating the horizontal position of
#' node labels. Must be one of \code{"left"}, \code{"center"}, or
#' \code{"right"}. Ignored if \code{label=FALSE}.
#' @param edge_width An integer representing the thickness of edges, or the
#' lines that connect nodes.
#' @param edge_color A valid hex color code as a string for the edge lines.
#' @template -browser
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE}, additional output
#' is printed to the console and the sequences, resulting clusters (convergence
#' groups), and XML file for Cytoscape are written to a folder in the working
#' directory.
#' @import htmlwidgets
#' @examples
#' data(tcr)  # Packaged data set, a data.frame from a CSV file
#' tcr_reduced <- tcr[1:50, ]
#' @export
convergence <- function(d, seqs_col=NULL, background_color="#FFFFFF",
                        node_shape="ELLIPSE", border_width=2,
                        border_color="#161616", node_color="#0B94B1",
                        node_size=30, labels=TRUE, label_vertical_pos="middle",
                        label_horizontal_pos="center", edge_width=3,
                        edge_color="#2D2D2D", browser=FALSE, verbose=FALSE) {
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
  xml_file <- cytoscape_xml(clusters, row_num=largest_cluster, labels, verbose,
                            verbose_dir)
  
  xml <- XML::xmlRoot(XML::xmlTreeParse(xml_file))
  xml_string <- XML::toString.XMLNode(xml)
  
  # Forward options and the XML string to convergence.js using 'x'
  x <- list(
    xml_string=xml_string,
    background_color=background_color,
    node_shape=node_shape,
    border_width=border_width,
    border_color=border_color,
    node_color=node_color,
    node_size=node_size,
    label_vertical_pos=label_vertical_pos,
    label_horizontal_pos=label_horizontal_pos,
    edge_width=edge_width,
    edge_color=edge_color
  )
  
  # Add the XML file as an HTML dependency so it can get loaded in the browser.
  # The Cytoscape visualization doesn't use this (it uses the xml_string above),
  # but it can be useful to debug by simply looking at this file in the browser.
  convergence_xml <- htmltools::htmlDependency(
    name = "convergence_xml",
    version = "1.0",
    src = c(file=dirname(xml_file)),
    attachment = list(xml=basename(xml_file))
  )
  
  # Create widget
  htmlwidgets::createWidget(
    name = "convergence",
    x,
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
