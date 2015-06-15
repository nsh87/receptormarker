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
    stop(seqs_col_err, call.=FALSE)
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


calculate_canvas_size <- function(xml_file, rings) {
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
  } else if (num_elements <= 30) {
    base_size <- 1300
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
  # Need to add extra room if there are outer rings (100 pixels per ring works
  # for up to 3 rings, check if it still does if using more rings)
  if (!is.null(rings)) {
    canvas_size <- canvas_size + length(rings) * 100
  }
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
  # Add the bar info and tooltip text to each element
  # <clade>
  #   <name>CATSREWGEAYEQYF</name>
  #   <branch_length>0/141994136</branch_length>
  #   <chart>
  #     <content>3</content>  # This makes the bar length
  #   </chart>
  #   <annotation>
  #     <desc>Expansions: 3</desc>  # This is the Unitip tooltip
  #   </annotation>
  # </clade>
  node_parents <- lapply(named_nodes, function(n) { XML::xmlParent(n) })
  node_parents <- lapply(named_nodes, function(n) {
    node_name <- XML::xmlValue(n)
    parent <- XML::xmlParent(n)
    
    # Subtract one from each seq's frequency so frequencies == 1 don't get bars
    bar_length <- seq_frequencies[names(seq_frequencies) == node_name] - 1
    desc_text <- paste0(c("Expansions:", bar_length), collapse=" ")
    
    chart_el <- XML::newXMLNode("chart")
    bar_el <- XML::newXMLNode("content", bar_length)
    annotation_el <- XML::newXMLNode("annotation")
    desc_el <- XML::newXMLNode("desc", desc_text)
    
    chart_el <- XML::addChildren(chart_el, kids=list(bar_el), append=FALSE)
    annotation_el <- XML::addChildren(annotation_el, kids=list(desc_el))
    parent <- XML::addChildren(parent, kids=list(chart_el,
                                                 annotation_el), append=FALSE)
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


validate_rings <- function(rings, d) {
  tryCatch({
    # Rings should be either NULL or a named character or numeric vector
    if (is.null(rings)) {
      return()
    } else if (!(class(rings) %in% c("character", "numeric"))) {
      err <- paste0(c("The argument 'rings' must be a named vector,",
                      "or NULL if no rings are desired"), collapse=" ")
      stop(err, call.=FALSE)
    } else if (length(names(rings)[names(rings) != ""]) != length(rings)) {
      stop("All elements in the argument 'rings' must be named", call.=FALSE)
    } else if (any(is.na(rings))) {
      # You can't have rings <- c(seqs=NA)
      stop("Values in the argument 'rings' cannot be NA", call.=FALSE)
    } else if (class(d) != "data.frame") {
      err <- paste0(c("Argument 'd' expected to be a data.frame when using",
                      "the argument 'rings'"), collapse=" ")
      stop(err, call.=FALSE)
    } else if (ncol(d) <= 1) {
      err <- paste0(c("The data.frame 'd' must have more than one column",
                      "in order to use the argument 'rings'"), collapse=" ")
      stop(err, call.=FALSE)
    } else if (!(all(names(rings) %in% names(d)))) {
      # Check that all the names in 'rings' are column headers of 'd', and if
      # not tell user which one is not in 'd'
      for (i in c(1:length(rings))) {
        if (!(names(rings)[i] %in% names(d))) {
          err <- paste0(c("Validation of the argument 'rings' failed: the ",
                          "column name '", names(rings)[i], "' does not ",
                          "exist in input data 'd'"), collapse="")
          stop(err, call.=FALSE)
        }
      }
    } else {
      return()
    }
  },
  error = function(e) {
    stop(e)
  }
  )
}


validate_true_false <- function(arg_list) {
  for (a in arg_list) {
    if (!(a %in% c(TRUE, FALSE))) {
      err <- paste0("The argument '", names(a), "' must be TRUE or FALSE",
                    collapse="")
      stop(err, call.=FALSE)
    }
  }
}


index_with_wrap <- function(index, v) {
  indx <- index %% length(v)
  if (indx == 0 ) { indx = length(v) }
  list(name=names(v)[indx], val=v[indx])
}


add_phylo_outer_rings <- function(xml_file, seqs, d_clean, seqs_col,
                                  rings, condensed) {
  colors <- c(green="#82A538", orange="#B1903C", blue="#626DBC", pink="#B5598F",
              turqoise="#0DA765", yellow="#E9E700", purple="#64059C",
              tangerine="#DE531C")
  # If condensed == TRUE, the sequences in xml_file are unique; in order to get
  # attributes for each ring, for each sequence you have to average each
  # attribute across all duplicate sequences, which means two algorithms
  # (one if condensed, one if not)
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  unique_seqs <- unique(seqs)
  # For each sequence...
  for (seq in unique_seqs) {
    search_str <- paste0("//ns:name[text()='", seq, "']")  # Grab all XML nodes
    seq_nodes <- XML::getNodeSet(root, search_str, namespaces=ns)
    # Get all rows of data where sequence == seq
    seq_data_rows <- d_clean[d_clean[, seqs_col] == seq, ]
    # Need to add rings data to the parent (clade) of each XML node:
    # <clade>
    #   <name>CATSREWGEAYEQYF</name>
    #   <branch_length>0.141994136</branch_length>
    #   <chart>
    #     <outergroup1>green</outergroup1>  # This colors first ring green
    #     <outergroup4>pink</outergroup4>  # This colors fourth ring pink
    #   </chart>
    # </clade>
    node_parents <- lapply(seq_nodes, XML::xmlParent)
    lapply(seq_along(node_parents), function(y, i) {
      parent <- y[[i]]  # Current parent node
      rings_to_add <- vector(mode="list", length=length(rings))
      if (condensed == FALSE) {
        # Compare 'rings' to the attributes of the current row of data
        data_row <- seq_data_rows[i, ]
        # Get values of current row's attributes in named list
        attrs <- c(data_row[names(rings)])
        # For each attribute...
        for (j in c(1:length(attrs))) {
          # If the value in the data matches the criteria in 'rings'...
          if (!is.na(attrs[[j]]) && attrs[[j]] == rings[j]) {
            # Add an XML node to color in the outer ring 
            outer_ring_name <- paste0(c("outergroup", j), collapse="")
            outer_ring_el <- XML::newXMLNode(outer_ring_name,
                                             index_with_wrap(j, colors)[["name"]])
            rings_to_add[[j]] <- outer_ring_el
          }
        }
        # No <chart> node if condensed == FALSE, so add it now
        chart_el <- XML::newXMLNode("chart")
        rings_to_add <- rings_to_add[!is.na(rings_to_add)]
        if (length(rings_to_add) > 0) {
          chart_el <- XML::addChildren(chart_el, kids=rings_to_add)
          parent <- XML::addChildren(parent, chart_el)
        }
      } else if (condensed == TRUE) {
        for (j in c(1:length(rings))) {
          # If >50% of the rows have an attribute that matches 'rings' add ring
          this_attr <- seq_data_rows[, names(rings)[j]]
          n_matches <- sum(this_attr == rings[j], na.rm=TRUE)
          if (n_matches / length(this_attr) >= 0.5) {
            outer_ring_name <- paste0(c("outergroup", j), collapse="")
            outer_ring_el <- XML::newXMLNode(outer_ring_name,
                                             index_with_wrap(j, colors)[["name"]])
            rings_to_add[[j]] <- outer_ring_el
          }
        }
        # The <chart> node exists if condensed == TRUE
        chart_el <- XML::xmlElementsByTagName(parent, "chart")[[1]]
        rings_to_add <- rings_to_add[!is.na(rings_to_add)]
        if (length(rings_to_add) > 0) {
          chart_el <- XML::addChildren(chart_el, kids=rings_to_add)
          parent <- XML::addChildren(parent, chart_el)
        }
      }
      parent
    }
    , y=node_parents
    )
  }
  # Add rings styles. If condensed == FALSE there won't be a <styles> or
  # <charts> section yet. Have to create them both. <charts> section dictates
  # how the rings will be drawn and what data type they are. <style> matches
  # up with the colors ('orange', 'blue', etc.)
  phylogeny_set <- XML::getNodeSet(root,
                                   "/ns:phyloxml/ns:phylogeny", namespaces=ns)
  phylogeny <- phylogeny_set[[1]]
  # Outer rings color styles
  outergroups <- vector(mode="list", length=length(rings))
  for (q in c(1:length(rings))) {
    outergroup_name <- paste0(c("outergroup", q), collapse="")
    outergroups[[q]] <- XML::newXMLNode(outergroup_name, attrs=c(disjointed=q,
                                                               type="binary"))
  }
  # List names are not preserved in lapply...have to use indices
  color_styles <- lapply(seq_along(colors), function(y, n, i) {
    XML::newXMLNode(n[[i]], attrs=c(fill=y[[i]], opacity="#0.95",
                                    stroke="#DDDDDD"))  
  }
  , y=colors, n=names(colors))
  if (condensed == FALSE) {
    render_el <- XML::newXMLNode("render")
    charts_el <- XML::newXMLNode("charts")
    charts_el <- XML::addChildren(charts_el, kids=outergroups)
    styles_el <- XML::newXMLNode("styles")
    styles_el <- XML::addChildren(styles_el, kids=color_styles)
    render_el <- XML::addChildren(render_el, charts_el, styles_el)
    phylogeny <- XML::addChildren(phylogeny, render_el)
  } else if (condensed == TRUE) {
    charts_el <- XML::getNodeSet(root,
                                "/ns:phyloxml/ns:phylogeny/ns:render/ns:charts",
                                namespaces=ns)[[1]]
    charts_el <- XML::addChildren(charts_el, kids=outergroups)
    styles_el <- XML::getNodeSet(root,
                                "/ns:phyloxml/ns:phylogeny/ns:render/ns:styles",
                                namespaces=ns)[[1]]
    styles_el <- XML::addChildren(styles_el, kids=color_styles)
  }
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}


validate_not_null <- function(arg_list) {
  for (i in c(1:length(arg_list))) {
    a <- arg_list[[i]]
    if (is.null(a)) {
      err <- paste0("The argument '", names(arg_list)[[i]],
                    "' cannot be NULL", collapse="")
      stop(err, call.=FALSE)
    }
  }
}


#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(d, seqs_col=NULL, condense=FALSE, rings=NULL,
                         canvas_size="auto", font_size="auto", scale=TRUE,
                         browser=FALSE, verbose=FALSE, fast=FALSE) {
  
  check_muscle(level="stop")
  biopy_existence <- check_bio_python(level="warn")
  
  # Validate function parameters
  validate_not_null(list(d=d, condense=condense, canvas_size=canvas_size,
                         font_size=font_size, scale=scale, browser=browser,
                         verbose=verbose, fast=fast))
  validate_canvas_size(canvas_size)
  seqs <- extract_sequences (d, seqs_col)
  validate_sequences(seqs)
  validate_rings(rings, d)
  validate_true_false(c(condense=condense, scale=scale, browser=browser,
                        verbose=verbose, fast=fast))
  
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
  
  # Condense phylo if condense == TRUE
  if (condense == TRUE) {
    # Send the non-unique sequences to be able to count how many there are
    xml_file <- add_bars_to_condensed_phyloxml(xml_file, clean[["seqs"]])
  }
  
  # Add outer-rings to phylo if there are any requested
  if (!is.null(rings)) {
    xml_file <- add_phylo_outer_rings(xml_file, seqs, clean[["d"]], seqs_col,
                                      rings, condense)
  }
  
  # Also write the phyloxml to the verbose folder if the user wants it
  if (verbose == TRUE && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  
  # Calculate canvas size based on number of nodes in phylo.xml
  if (canvas_size == "auto") {
    canvas_size <- calculate_canvas_size(xml_file, rings) 
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
