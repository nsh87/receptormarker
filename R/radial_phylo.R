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


extract_sequences <- function(df, seqs_col) {
  seqs_col_err <- "The argument 'df' and/or 'seqs_col' is invalid" 
  tryCatch({
    if (length(seqs_col) != 1) {
      stop(seqs_col_err, call.=FALSE)
    }
    if (typeof(seqs_col) == "character" && seqs_col %in% names(df)) {
      seqs <- as.character(df[, seqs_col])
    } else if (seqs_col == floor(seqs_col) && length(seqs_col) == 1) {
      seqs <- as.character(df[, seqs_col])
    } else {
      stop(seqs_col_err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(seqs_col_err, call.=FALSE)
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


clean_df <- function(df, seqs_col, verbose, verbose_dir) {
  df_clean <- df[df[, seqs_col] != "", ]  # Remove rows with no sequences
  df_clean <- df_clean[complete.cases(df_clean[, seqs_col]), ]  # Remove NA rows
  # Write the sequences if the user wants them
  if (verbose == TRUE) {
    seqs_fasta <- with(df_clean,
                       paste0(">", df_clean[, seqs_col], "\n",
                              df_clean[, seqs_col],
                              collapse="\n"
                              )
    )
    seqs_file <- tempfile(pattern="sequences-", tmpdir=verbose_dir,
                          fileext=".txt")
    write(seqs_fasta, seqs_file)
  }
  df_clean
}


msa <- function(seqs, verbose, verbose_dir) {
  seqs_biostring <- Biostrings::AAStringSet(seqs)
  names(seqs_biostring) <- seqs
  if (verbose == TRUE) message("MUSCLE multiple sequence alignment:")
  ms_alignment <- muscle::muscle(stringset=seqs_biostring, quiet=!verbose)
  # Write the alignment if the user wants it
  if (verbose == TRUE) {
    aa_str_set <- as(ms_alignment, "AAStringSet")
    msa_file <- tempfile(pattern="msa-", tmpdir=verbose_dir, fileext=".fasta")
    Biostrings::writeXStringSet(aa_str_set, file=msa_file)
  }
  ms_alignment
}


compute_dist_matrix <- function(ms_alignment, seqs) {
  # Need to turn the MSA into class 'alignment'
  alignment <- seqinr::as.alignment(nb=nrow(ms_alignment),
                            nam=seqs,
                            seq=as.character(ms_alignment, use.names=FALSE)) 
  dist_matrix <- as.matrix(seqinr::dist.alignment(x=alignment,
                                                  matrix="identity"))
}


create_tree <- function(dist_matrix, verbose, verbose_dir) {
  dist_tree <- ape::bionj(dist_matrix)
  phylo_tree <- ape::as.phylo(dist_tree)
  newick_file <- tempfile(pattern="tree-", tmpdir=tempdir(), fileext=".newick")
  ape::write.tree(phy=phylo_tree, file=newick_file)
  # Copy the newick file from the tmp dir to verbose dir if the user wants it
  if (verbose == TRUE && file.exists(newick_file)) {
    file.copy(newick_file, verbose_dir, copy.date=TRUE)
  }
  newick_file
}


newick_to_phyloxml <- function(newick_file, verbose, verbose_dir) {
  # Create tempdir for phyloxml file
  phyloxml_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(phyloxml_tmpdir)  # htmlwidgets copies entire dir to browser
  
  xml_file <- tempfile(pattern="phyloxml-", tmpdir=phyloxml_tmpdir,
                       fileext=".xml")
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
  # Also write it to the verbose folder if the user wants it
  if (verbose == TRUE && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir, copy.date=TRUE)
  }
  xml_file
}


calculate_canvas_size <- function(xml_file) {
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  num_elements <- length(XML::xpathApply(root, "//ns:name", namespaces=ns))
  # Similarly: xpathApply(doc, "/ns:phyloxml//ns:name", xmlValue, namespaces=ns)
  if (num_elements == 0) {
    stop("Cannot generate phylogram: no sequences to plot", call.=FALSE)
  } else if (num_elements <= 50) {
    return(900)
  } else if (num_elements <= 185) {
    return(1300)
  } else if (num_elements <= 1000) {
    return(num_elements * 7)
  } else {
    performance_warning <- paste0("Performance of the phylogram plot might",
                                  "begin to degrade with >1000 sequences",
                                  collapse=" ")
    message(performance_warning)
    return(num_elements * 7)
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
radial_phylo <- function(df, seqs_col, canvas_size="auto", font_size="auto",
                         scale=TRUE, browser=FALSE, verbose=FALSE) {
  
  # Validate function parameters
  validate_canvas_size(canvas_size)
  seqs <- extract_sequences (df, seqs_col)
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
  df_clean <- clean_df(df, seqs_col, verbose, verbose_dir)
  seqs <- as.character(df_clean[, seqs_col])
  # Step 2: Do a multiple sequence alignment (MSA)
  ms_alignment <- msa(seqs, verbose, verbose_dir)
  # Step 3: Compute pairwise distances of the aligned sequences
  dist_matrix <- compute_dist_matrix(ms_alignment, seqs)
  # Step 4: Calculate a distance tree and write it as .newick
  newick_file <- create_tree(dist_matrix, verbose, verbose_dir)
  # Step 5: Convert the .newick to phylo.xml
  xml_file <- newick_to_phyloxml(newick_file, verbose, verbose_dir)
  
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
