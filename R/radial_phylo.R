validate_canvas_size <- function(canvas_size) {
  err <- "The argument 'canvas_size' is invalid"
  canvas_size_options <- "auto"
  tryCatch({
    if (length(canvas_size) != 1) {
      stop(err, call.=FALSE)
    }
    if (!is.element(canvas_size, canvas_size_options) &&
        canvas_size != floor(canvas_size)) {
      stop(err, call.=FALSE)
    } else if (canvas_size < 1) {
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(err, call.=FALSE)
  }
  )
}


extract_sequences <- function(d, seqs_col) {
  # Validation checks
  if (is.null(d)) {
    stop("Data cannot be NULL", call.=FALSE)
  } else if (!(class(d) %in% c("data.frame", "character"))) {
    stop("Data must be a data.frame or a character vector", call.=FALSE)
  } else if (class(d) == "data.frame" &&
             !(class(seqs_col) %in% c("character", "numeric"))) {
    err <- paste0(c("Argument 'seqs_col' must be a column index or a column",
                    "header when data is a data.frame, or NULL if data is",
                    "a vector"),
                  collapse = " ")
    stop(err, call.=FALSE)
  } else if (class(d) == "character" && !is.null(seqs_col)) {
    err <- "Argument 'seqs_col' expected to be NULL when data is a vector"
    stop(err, call.=FALSE)
  } else if (class(d) == "data.frame" && length(seqs_col) != 1) {
    err <- "Argument 'seqs_col' must represent a single column, i.e. length = 1"
    stop(err, call.=FALSE)
  } else if (class(d) == "data.frame" &&
             class(seqs_col) == "numeric" &&
             seqs_col > ncol(d)) {
    stop("Argument 'seqs_col' is greater than the number of columns of data",
         call.=FALSE)
  } else if (typeof(seqs_col) == "character" && !(seqs_col %in% names(d))) {
    err <- paste0(c("Data does not have column named '", seqs_col, "'"),
                  collapse="")
    stop(err, call.=FALSE)
  }
  
  # Return the sequences
  tryCatch({
    # If d is single-column data frame return the sequences
    if (class(d) == "data.frame" && dim(d)[2] == 1) {
      seqs <- as.character(d[, 1])
      return(seqs)
    # Otherwise check if d is a vector, return it if so
    } else if (class(d) == "character") {
      return(d)
    # If data is >1D data.frame return the sequences column by header name
    } else if (typeof(seqs_col) == "character" && seqs_col %in% names(d)) {
      seqs <- as.character(d[, seqs_col])
      return(seqs)
    # Or if data is >1D data.frame return the sequences column by index
    } else if (seqs_col == floor(seqs_col) && length(seqs_col) == 1) {
      seqs <- as.character(d[, seqs_col])
      return(seqs)
    } else {
      err <- paste0(c("Unexpected error retrieving sequences. Please check",
                      "paraameters 'd' and 'seqs_col'."), collapse=" ")
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(e)
  }
  )
}


validate_sequences <- function(seqs) {
  # Make sure sequences are only alpha characters
  seqs_col_err <- "Sequences must only contain characters from A-Z"
  g <- grepl("[^A-Za-z]", as.character(seqs))
  if (sum(g) > 0) {
    stop(seqs_col_err, cal.=FALSE)
  }
}


clean_d <- function(d, seqs_col, verbose, verbose_dir) {
  if (class(d) == "data.frame" && dim(d)[2] > 1) {
    d_clean <- d[d[, seqs_col] != "", ]  # Remove rows with no sequences
    d_clean <- d_clean[complete.cases(d_clean[, seqs_col]), ]  # Remove NA rows
    d_clean <- d_clean[with(d_clean, order(d_clean[, seqs_col])), ]
    seqs <- as.character(d_clean[, seqs_col])
  } else if (class(d) == "data.frame" && dim(d)[2] == 1) {
    d_clean <- d[d != ""]  # Remove blank sequences
    d_clean <- d_clean[!is.na(d_clean)]  # Remove NA's
    d_clean <- sort(d_clean)
    seqs <- as.character(d_clean)
  } else if (class(d) == "character") {
    d_clean <- d[d != ""]  # Remove blank sequences
    d_clean <- d_clean[!is.na(d_clean)]  # Remove NA's
    d_clean <- sort(d_clean)
    seqs <- as.character(d_clean)
  }
  # Write the sequences if the user wants them
  if (verbose == TRUE) {
    if (class(d) == "data.frame") {
      seqs_fasta <- with(d_clean,
                         paste0(">", seqs, "\n",
                                seqs,
                                collapse="\n"
                                )
      )
    } else if (class(d) == "character") {
      seqs_fasta <- paste0(">", seqs, "\n",
                           seqs,
                           collapse="\n")
    }
    seqs_file <- tempfile(pattern="sequences-", tmpdir=verbose_dir,
                          fileext=".txt")
    write(seqs_fasta, seqs_file)
  }
  list("seqs"=seqs, "d"=d_clean)
}


msa <- function(seqs, dedupe_hash, verbose, verbose_dir) {
  seqs_biostring <- Biostrings::AAStringSet(seqs)
  # Create unique names for the sequences to be written to file so Biopython
  # doesn't complain about duplicate sequence names. We use the counts of each
  # sequence to add a suffix (e.g. if 'CASHT' if found 3 times they will become
  # 'CASHT.1', 'CASHT.2', 'CASHT.3'); and we prepend that with a hash so that
  # we can reliably find our own dedupe sequence and remove it later. If we
  # didn't have the has and just tried to remove '.1', '.2', and '.3' later we
  # might end up removing part of the user's sequence name.
  seqs_rle <- rle(seqs)
  seqs_unique <- paste0(rep(seqs_rle[['values']], times=seqs_rle[['lengths']]),
                        ".", dedupe_hash, "-",
                        unlist(lapply(seqs_rle[['lengths']], seq_len)))
  names(seqs_biostring) <- seqs_unique
  if (verbose == TRUE) print("MUSCLE multiple sequence alignment:")
  ms_alignment <- muscle::muscle(stringset=seqs_biostring, quiet=!verbose)
  # Write the alignment to file in case Biopython needs it
  aa_str_set <- as(ms_alignment, "AAStringSet")
  msa_file <- tempfile(pattern="msa-", tmpdir=tempdir(), fileext=".fasta")
  Biostrings::writeXStringSet(aa_str_set, file=msa_file)
  # Copy the alignment from the tmp dir to verbose dir if the user wants it
  if (verbose == TRUE && file.exists(msa_file)) {
    file.copy(msa_file, verbose_dir)
  }
  list(as_string=ms_alignment, file=msa_file)
}


compute_dist_matrix <- function(ms_alignment, seqs) {
  # Need to turn the MSA into class 'alignment'
  alignment <- seqinr::as.alignment(nb=nrow(ms_alignment),
                            nam=seqs,
                            seq=as.character(ms_alignment, use.names=FALSE)) 
  dist_matrix <- as.matrix(seqinr::dist.alignment(x=alignment,
                                                  matrix="identity"))
}


create_tree <- function(dist_matrix) {
  dist_tree <- ape::bionj(dist_matrix)
  phylo_tree <- ape::as.phylo(dist_tree)
  newick_file <- tempfile(pattern="tree-", tmpdir=tempdir(), fileext=".newick")
  ape::write.tree(phy=phylo_tree, file=newick_file)
  newick_file
}


phyloxml_temp <- function() {
  # Create tempdir for phyloxml file
  phyloxml_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(phyloxml_tmpdir)  # htmlwidgets copies entire dir to browser
  
  xml_file <- tempfile(pattern="phyloxml-", tmpdir=phyloxml_tmpdir,
                       fileext=".xml")
}


newick_to_phyloxml <- function(newick_file, verbose, verbose_dir) {
  xml_file <- phyloxml_temp()
  forester <- system.file("java", "forester_1038.jar", package="receptormarker")
  system(sprintf(
    "java -cp %s org.forester.application.phyloxml_converter -f=nn -ni %s %s",
    forester,
    newick_file,
    xml_file
    ),
    ignore.stdout=!verbose,
    ignore.stderr=!verbose
  )
  # Also write the phyloxml to the verbose folder if the user wants it
  if (verbose == TRUE && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  xml_file
}


phyloxml_via_biopython <- function(ms_alignment, verbose) {
  xml_file <- phyloxml_temp()
  phyloxml_from_msa <- system.file("py", "phyloxml_from_msa.py",
                                   package="receptormarker")
  system(sprintf(
    "python %s --msa %s --dest %s",
    phyloxml_from_msa,
    ms_alignment,
    xml_file
    ),
    ignore.stdout=!verbose,
    ignore.stderr=!verbose
  )
  xml_file
}


remove_phyloxml_hash <- function(xml_file, hash) {
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  named_nodes <- XML::getNodeSet(root, "//ns:name", namespaces=ns)
  # Grab the separated period + hash + rest of chars until end of str
  regex_find <- paste0("\\.", hash, ".*$")  
  lapply(named_nodes, function(n) {
    XML::xmlValue(n) <- gsub(regex_find, "", XML::xmlValue(n))
  })
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}


calculate_canvas_size <- function(xml_file) {
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  named_nodes <- XML::xpathApply(root, "//ns:name", XML::xmlValue,
                                 namespaces=ns)
  named_nodes <- as.character(named_nodes)
  # Biopython adds a name "Inner123" to nodes that don't get shown. Need to
  # not count these when determining how many names are on phylogram
  named_nodes <- named_nodes[grepl("^[^Inner]", named_nodes)]
  num_elements <- length(named_nodes)
  # Similarly: xpathApply(doc, "/ns:phyloxml//ns:name", xmlValue, namespaces=ns)
  if (num_elements == 0) {
    stop("Cannot generate phylogram: no sequences to plot", call.=FALSE)
  } else if (num_elements <= 50) {
    base_size <- 1100
  } else if (num_elements <= 185) {
    base_size <- 1300
  } else if (num_elements <= 1000) {
    base_size <- num_elements * 7
  } else {
    performance_warning <- paste0(c("Depending on your computer, performance",
                                    "might begin to degrade when plotting",
                                    ">1000 sequences"),
                                  collapse=" ")
    warning(performance_warning, call.=FALSE)
    base_size <- num_elements * 7
  }
  # Need to take into account length of names on the phylogram plot
  longest_sequence <- max(nchar(named_nodes))
  # Your base_size should allow for up to 20 char sequence names; beyond 15
  # chars you need to expand the canvas to give them more room.
  quotient <- longest_sequence %/% 20
  remainder <- 0
  if (quotient > 0) {
    remainder <- longest_sequence %% 20
  }
  extra_room <- quotient * 500 + remainder * 60
  canvas_size <- base_size + extra_room
  canvas_size
}


add_bars_to_condensed_phyloxml <- function(xml_file, seqs) {
  seq_frequencies <- table(seqs)
  
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  # Get all elements that don't start with Biopython's "InnerXXX"
  named_nodes <- XML::getNodeSet(root,
                                 "//ns:name[not(starts-with(text(), 'Inner'))]",
                                 namespaces=ns)
  # Add the bar info to each element
  node_parents <- lapply(named_nodes, function(n) { XML::xmlParent(n) })
  node_parents <- lapply(named_nodes, function(n) {
    node_name <- XML::xmlValue(n)
    chart_el <- XML::newXMLNode("chart")
    # Subtract one from each seq's frequency so frequencies == 1 don't get bars
    bar_length <- seq_frequencies[names(seq_frequencies) == node_name] - 1
    bar_el <- XML::newXMLNode("content", bar_length)
    chart_el <- XML::addChildren(chart_el, kids=list(bar_el), append=FALSE)
    parent <- XML::xmlParent(n)
    parent <- XML::addChildren(parent, kids=list(chart_el), append=FALSE)
    parent
  })
  # Add the bar styles, which need to be appended to the phylogeny
  phylogeny_set <- XML::getNodeSet(root,
                                   "/ns:phyloxml/ns:phylogeny", namespaces=ns)
  phylogeny <- phylogeny_set[[1]]
  # <render>
  #  <charts>
  #    <content fill="#666" type="bar" width="0.2"/>
  #  </charts>
  #  <styles>
  #    <barChart fill="#333" stroke-width="0"/>
  #  </styles>
  # </render>
  render_el <- XML::newXMLNode("render")
  charts_el <- XML::newXMLNode("charts")
  content_el <- XML::newXMLNode("content", attrs=c(fill="#666",
                                                   type="bar",
                                                   width="0.2"))
  charts_el <- XML::addChildren(charts_el, kids=list(content_el))
  styles_el <- XML::newXMLNode("styles")
  barchart_el <- XML::newXMLNode("barChart", attrs=c(fill="#333",
                                                     "stroke-width"="0"))
  styles_el <- XML::addChildren(styles_el, kids=list(barchart_el))
  render_el <- XML::addChildren(render_el, kids=list(charts_el, styles_el))
  phylogeny <- XML::addChildren(phylogeny, kids=list(render_el))
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}


#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(d, seqs_col=NULL, condense=FALSE, canvas_size="auto",
                         font_size="auto", scale=TRUE, browser=FALSE,
                         verbose=FALSE, fast=FALSE) {
  
  check_muscle(level="stop")
  biopy_existence <- check_bio_python(level="warn")
  
  # Validate function parameters
  validate_canvas_size(canvas_size)
  seqs <- extract_sequences (d, seqs_col)
  validate_sequences(seqs)
  
  # Not necessary to have as func parameters; these will get set automatically
  width <- NULL
  height <- NULL
  
  # Create verbose dir
  if (verbose == TRUE) {
    verbose_dir <- tempfile("radial_phylo-", tmpdir=getwd(), fileext="")
    dir.create(verbose_dir)
  }
  
  # Step 1: Clean the data.frame and get the cleaned sequences
  clean <- clean_d(d, seqs_col, verbose, verbose_dir)
  seqs <- clean[["seqs"]]
  if (condense == TRUE) {
    seqs <- unique(seqs)
  }
  # Step 2: Do a multiple sequence alignment (MSA)
  dedupe_hash <- paste0(sample(c(0:9, letters), 10, replace=TRUE), collapse="")
  ms_alignment <- msa(seqs, dedupe_hash, verbose, verbose_dir)
  if (is.null(biopy_existence) || fast == TRUE) {
    # Step 3: Compute pairwise distances of the aligned sequences
    dist_matrix <- compute_dist_matrix(ms_alignment[["as_string"]], seqs)
    # Step 4: Calculate a distance tree and write it as .newick
    newick_file <- create_tree(dist_matrix)
    # Step 5: Convert the .newick to phylo.xml
    xml_file <- newick_to_phyloxml(newick_file, verbose, verbose_dir)
  } else {
    xml_file_with_hash <- phyloxml_via_biopython(ms_alignment[["file"]],
                                                 verbose)
    xml_file <- remove_phyloxml_hash(xml_file_with_hash, dedupe_hash)
  }
  
  # Add annotations to phylo.xml
  if (condense == TRUE) {
    # Send the non-unique sequences to be able to count how many there are
    xml_file <- add_bars_to_condensed_phyloxml(xml_file, clean[["seqs"]])
  }
  
  # Also write the phyloxml to the verbose folder if the user wants it
  if (verbose == TRUE && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  
  # Calculate canvas size based on number of nodes in phylo.xml
  if (canvas_size == "auto") {
    canvas_size <- calculate_canvas_size(xml_file) 
  }
  
  # Forward options to radial_phylo.js using 'x'
  x <- list(
    canvas_size = canvas_size,
    scale = scale
  )
  
  # Add the phyloxml as an HTML dependency so it can get loaded in the browser
   phyloxml <- htmltools::htmlDependency(
     name = "phyloxml",
     version = "1.0",
     src = c(file=dirname(xml_file)),
     attachment = list(xml=basename(xml_file))
   )
   
  # Create widget
  htmlwidgets::createWidget(
    name = "radial_phylo",
    x,
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(
      padding = 22,
      viewer.suppress = browser,
      browser.fill = TRUE
    ),
    package = "receptormarker",
    dependencies = phyloxml
  )
  
}


#' Widget output function for use in Shiny
#'
#' @export
# nolint start
radial_phyloOutput <- function(outputId, width = "100%", height = "400px"){
  shinyWidgetOutput(outputId, "radial_phylo", width, height,
                    package="receptormarker")
# nolint end
}


#' Widget render function for use in Shiny
#'
#' @export
# nolint start
renderRadial_phylo <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr) # force quoted
  shinyRenderWidget(expr, radial_phyloOutput, env, quoted = TRUE)
# nolint end
}
