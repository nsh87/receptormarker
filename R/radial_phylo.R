#' @title Validate a font size to be used with a radial phylogram
#' @description An internal function that raises an error if a font size is not
#' an integer between 1 and 99.
#' @param font_size A font size
#' @keywords internal
validate_font_size <- function(font_size) {
  if (any(is.na(font_size))) {
    stop("Argument 'font_size' must be a real value", call.=FALSE)
  } else if (class(font_size) != "numeric") {
    stop("Argument 'font_size' must be an integer", call.=FALSE)
  } else if (length(font_size) != 1) {
    stop("Argument 'font_size' must be a single integer", call.=FALSE)
  } else if(font_size == 0) {
    stop("Argument 'font_size' must be greater than 0", call.=FALSE)
  } else if (length(grep("^[0-9]{1,3}$", font_size)) != 1 ||
             grep("^[0-9]{1,2}$", font_size) != 1) {
    stop("Argument 'font_size' must be an integer between 1 and 99",
         call.=FALSE)
  }
}


#' @title Validate the tree margin to be used with a radial phylogram
#' @description An internal function that raises an error if a tree margin
#' is not a number between 0 and 1.
#' @param tree_margin A tree margin
#' @keywords internal
validate_tree_margin <- function(tree_margin) {
  if (any(is.na(tree_margin))) {
    stop("Argument 'tree_margin' must be a real value", call.=FALSE)
  } else if (class(tree_margin) != "numeric") {
    stop("Argument 'tree_margin' must be a number", call.=FALSE)
  } else if (length(tree_margin) != 1) {
    stop("Argument 'tree_margin' must be a single number", call.=FALSE)
  } else if (!(tree_margin > 0)) {
    stop("Argument 'tree_margin' must be greater than 0 and less than 1",
         call.=FALSE)
  } else if (!(tree_margin < 1)) {
    stop("Argument 'tree_margin' must be greater than 0 and less than 1",
         call.=FALSE)
  }
}


#' @title Validate a ring thickness to be used with a radial phylogram
#' @description An internal function that raises an error if a ring
#' thickness is not an integer between 1 and 99.
#' @param ring_thickness A ring thickness
#' @keywords internal
validate_ring_thickness <- function(ring_thickness) {
  if (any(is.na(ring_thickness))) {
    stop("Argument 'ring_thickness' must be a real value", call.=FALSE)
  } else if (class(ring_thickness) != "numeric") {
    stop("Argument 'ring_thickness' must be an integer", call.=FALSE)
  } else if (length(ring_thickness) != 1) {
    stop("Argument 'ring_thickness' must be a single integer", call.=FALSE)
  } else if(ring_thickness == 0) {
    stop("Argument 'ring_thickness' must be greater than 0", call.=FALSE)
  } else if (length(grep("^[0-9]{1,3}$", ring_thickness)) != 1 ||
             grep("^[0-9]{1,2}$", ring_thickness) != 1) {
    stop("Argument 'ring_thickness' must be an integer between 1 and 99",
         call.=FALSE)
  }
}


#' @title Validate arc size to be used with a radial phylogram
#' @description An internal function that raises an error if an arc size is not
#' an integer between 1 and 359.
#' @param arc An arc size in degrees.
#' @keywords internal
validate_arc <- function(arc) {
  if (any(is.na(arc))) {
    stop("Argument 'arc' must be a real value", call.=FALSE)
  } else if (!(class(arc) %in% c("character", "numeric"))) {
    stop("Argument 'arc' must be an integer or 'auto'", call.=FALSE)
  } else if (length(arc) != 1) {
    stop("Argument 'arc' must be a single integer or 'auto'", call.=FALSE)
  } else if (class(arc) == "character" && arc != "auto") {
    stop("Argument 'arc' must be a single integer or 'auto'", call.=FALSE)
  } else if (arc == "auto") {
    return()
  } else if(arc == 0) {
    stop("Argument 'arc' must be greater than 0", call.=FALSE)
  } else if (length(grep("^[0-9]{1,3}$", arc)) != 1) {
    stop("Argument 'arc' must be an integer between 1 and 359", call.=FALSE)
  } else if (!(arc > 0 && arc < 360)) {
    stop("Argument 'arc' must be an integer between 1 and 359", call.=FALSE)
  }
}


#' @title Validate the distance calculator to use on sequence alignment
#' @description The distance calculator that BioPython creates uses any one of a
#' number of distance matrices. This function validates that the user input
#' a single one of these matrices to use for distance calculation.
#' @param dist A string representing the distance matrix to use for distance
#' calculations.
#' @keywords internal
validate_distance <- function(dist) {
  distances <- c("benner6", "benner22", "benner74", "blosum100",
                 "blosum30", "blosum35", "blosum40", "blosum45", "blosum50",
                 "blosum55", "blosum60", "blosum62", "blosum65", "blosum70",
                 "blosum75", "blosum80", "blosum85", "blosum90", "blosum95",
                 "feng", "fitch", "genetic", "gonnet", "grant", "ident",
                 "johnson", "levin", "mclach", "miyata", "nwsgappep", "pam120",
                 "pam180", "pam250", "pam30", "pam300", "pam60", "pam90", "rao",
                 "risler", "structure")
  if (any(is.na(dist))) {
    stop("Argument 'dist' must be a real value", call.=FALSE)
  } else if (class(dist) != "character") {
    stop("Argument 'dist' must be a string", call.=FALSE)
  } else if (length(dist) != 1) {
    stop("Argument 'dist' must be a single string", call.=FALSE)
  } else if (!(dist %in% distances)) {
    stop("Argument 'dist' is not a valid distance matrix", call.=FALSE) 
  }
}


#' @title Validate a radial phylogram's canvas size
#' @description An internal function that creates an error if \code{canvas_size}
#' is invalid.
#' @param canvas_size A \code{canvas_size} from the \code{\link{radial_phylo}}
#' function.
#' @keywords internal
validate_canvas_size <- function(canvas_size) {
  err <- "The argument 'canvas_size' is invalid"
  canvas_size_options <- "auto"
  tryCatch({
    if (length(canvas_size) != 1) {
      stop(err, call.=FALSE)
    } else if (!(class(canvas_size) %in% c("character", "numeric"))) {
      stop(err, call.=FALSE)
    } else if (!is.element(canvas_size, canvas_size_options) &&
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


#' @title Validate the data, sequence column, and sequences for a phylogram
#' @description An internal function that raises an error if the arguments are
#' invalid, the data cannot be subset using the provided sequence column, or the
#' sequences extracted are invalid.
#' @template -d
#' @template -seqs_col
#' @keywords internal
validate_d_seqs <- function(d, seqs_col) {
  tryCatch({
    # First check the arguments:
    if (!(class(d) %in% c("data.frame", "character"))) {
      err <- "Argument 'd' must be a data.frame or a character vector"
      stop(err, call.=FALSE)
    } else if (class(d) == "data.frame" &&
               !(class(seqs_col) %in% c("character", "numeric"))) {
      err <- paste0(c("Argument 'seqs_col' must be a column index or a column",
                      "name when 'd' is a data.frame, or NULL if 'd' is",
                      "a vector"),
                    collapse = " ")
      stop(err, call.=FALSE)
    } else if (class(d) == "character" && !is.null(seqs_col)) {
      err <- "Argument 'seqs_col' expected to be NULL when data is a vector"
      stop(err, call.=FALSE)
    } else if (class(d) == "data.frame" && length(seqs_col) != 1) {
      err <- paste0(c("Argument 'seqs_col' must represent a single",
                      "column, i.e. length = 1"), collapse=" ")
      stop(err, call.=FALSE)
    } else if (class(d) == "data.frame" && class(seqs_col) == "numeric" &&
               seqs_col > ncol(d)) {
      stop("Argument 'seqs_col' is greater than the number of columns of data",
           call.=FALSE)
    } else if (typeof(seqs_col) == "character" && !(seqs_col %in% names(d))) {
      err <- paste0(c("Data does not have column named '", seqs_col, "'"),
                    collapse="")
      stop(err, call.=FALSE)
      # Now check the sequences:
    } else if (class(d) == "data.frame" && dim(d)[2] == 1) {
      # If d is single-column data frame...
      seqs <- as.character(d[, 1])
      validate_sequences(seqs)
      return()
    # If d is a vector...
    } else if (class(d) == "character") {
      validate_sequences(d)
      return()
    # If dis >1D data frame...
    } else if (typeof(seqs_col) == "character" && seqs_col %in% names(d)) {
      seqs <- as.character(d[, seqs_col])
      validate_sequences(seqs)
      return()
    # If d is >1D data frame...
    } else if (seqs_col == floor(seqs_col) && length(seqs_col) == 1) {
      seqs <- as.character(d[, seqs_col])
      validate_sequences(seqs)
      return()
    } else {
      err <- paste0(c("Unexpected error retrieving sequences. Please check ",
                      "parameters 'd' and 'seqs_col'."), collapse="")
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(e)
  }
  )
}


#' @title Validate the \code{color_wheel} for a phylogram
#' @description An internal function that raises an error if \code{color_wheel}
#' is not a named character vectors or if the hex codes are not valid.
#' @param color_wheel The \code{color_wheel} argument supplied to make the
#' radial phylogram.
#' @keywords internal
validate_color_wheel <- function(color_wheel) {
  tryCatch({
    if (class(color_wheel) != "character") {
      err <- paste0(c("The argument 'color_wheel' must be a named character ",
                      "vector"), collapse="")
      stop(err, call.=FALSE)
    } else if (length(
        names(color_wheel)[names(color_wheel) != ""]) != length(color_wheel)) {
      err <- paste0(c("All elements in the argument 'color_wheel' must be ",
                      "named", collapse=""))
      stop(err, call.=FALSE)
    } else if (any(is.na(color_wheel))) {
      # You can't have rings <- c(black=NA)
      stop("Values in the argument 'color_wheel' cannot be NA", call.=FALSE)
    } else if (any(grepl("^#[0-9a-zA-Z]{6}$", color_wheel)) == FALSE) {
      err <- paste0(c("Values in the argument 'color_wheel' must be valid hex ",
                      "color codes", collapse=""))
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(e)
  }
  )
}


#' @title Validate the user-supplied argument \code{rings} for a phylogram
#' @description An internal function that raises an error if \code{rings} is not
#' a named character vector, if its names do not appear in \code{d} as column
#' names, if \code{d} is not a data frame, and more.
#' @param rings A user-supplied named character vector representing rings to add
#' to a phylogram.
#' @template -d
#' @keywords internal
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
    } else if ("all" %in% rings && length(rings) > 1) {
        err <- paste0(c("Only a single column can be annotated with rings ",
                        "when using the value 'all'"), collapse="")
        stop(err, call.=FALSE)
    } else {
      return()
    }
  },
  error = function(e) {
    stop(e)
  }
  )
}


#' @title Perform multiple sequence alignment on sequences
#' @description An internal function that utilizes Bioconductor's MUSCLE for
#' MSA. Writes the alignment in addition to returning it so the file can be
#' loaded by Biopython.
#' @template -seqs
#' @param verbose \code{TRUE} or \code{FALSE}, depending on whether or not the
#' output of the MSA should be printed and the alignment file should be copied
#' to the \code{verbose_dir}.
#' @template -verbose_dir
#' @param dedupe_hash A random string to use to deduplicate sequence names when
#' writing the MSA to FASTA. This is necessary to prevent Biopython from
#' complaining about duplicate names when reading the MSA. By passing this value
#' in, it can later be removed from any subsequent files, such as the final
#' PhyloXML.
#' @return A named list containing the multiple sequence alignment
#' \code{"as_string"} (a \code{MultipleAlignment} class) and a path to the MSA
#' FASTA \code{"file"}.
#' @keywords internal
msa <- function(seqs, verbose, verbose_dir, dedupe_hash=NULL) {
  seqs_biostring <- Biostrings::AAStringSet(seqs)
  seqs_rle <- rle(seqs)
  if (!is.null(dedupe_hash)) {
      # Use the counts of each sequence to add a suffix (e.g. if 'CASHT' if
      # found 3 times they will become 'CASHT.1', 'CASHT.2', 'CASHT.3'); and we
      # prepend that with the dedupe_hash so that we can reliably find our own
      # dedupe sequence and remove it later. If we didn't have the hash and just
      # tried to remove '.1', '.2', and '.3' later we might end up removing part
      # of the user's sequence name.
      seqs_unique <- paste0(rep(seqs_rle[["values"]],
                                times=seqs_rle[["lengths"]]),
                            ".", dedupe_hash, "-",
                            unlist(lapply(seqs_rle[["lengths"]], seq_len)))
  } else seqs_unique <- seqs
  names(seqs_biostring) <- seqs_unique
  if (verbose) message("MUSCLE multiple sequence alignment:")
  ms_alignment <- muscle::muscle(stringset=seqs_biostring, quiet=!verbose)
  # Write the alignment to file in case Biopython needs it
  aa_str_set <- as(ms_alignment, "AAStringSet")
  msa_file <- tempfile(pattern="msa-", tmpdir=tempdir(), fileext=".fasta")
  Biostrings::writeXStringSet(aa_str_set, file=msa_file)
  # Copy the alignment from the tmp dir to verbose dir if the user wants it
  if (verbose && file.exists(msa_file)) {
    file.copy(msa_file, verbose_dir)
  }
  list(as_string=ms_alignment, file=msa_file)
}


#' @title Compute a distance matrix from a MSA
#' @description An internal function that uses the \emph{seqinr} package to
#' compute an identity distance matrix.
#' @template -ms_alignment
#' @template -seqs
#' @return A distance matrix of class \code{dist}.
#' @keywords internal
compute_dist_matrix <- function(ms_alignment, seqs) {
  # Need to turn the MSA into class 'alignment'
  alignment <- seqinr::as.alignment(nb=nrow(ms_alignment), nam=seqs,
                            seq=as.character(ms_alignment, use.names=FALSE)) 
  dist_matrix <- as.matrix(seqinr::dist.alignment(x=alignment,
                                                  matrix="identity"))
}


#' @title Create a phylogenetic tree
#' @description An internal function that uses the \emph{ape} package and the
#' BIONJ method (an improved NJ, or neighbor joining, method) to create a
#' phylogenetic tree.
#' @param dist_matrix A distance matrix of class \code{dist}, likely created by
#' the \code{\link{compute_dist_matrix}} function.
#' @return A path to a Newick file representing the phylogenetic tree.
#' @keywords internal
create_tree <- function(dist_matrix) {
  dist_tree <- ape::bionj(dist_matrix)
  phylo_tree <- ape::as.phylo(dist_tree)
  newick_file <- tempfile(pattern="tree-", tmpdir=tempdir(), fileext=".newick")
  ape::write.tree(phy=phylo_tree, file=newick_file)
  newick_file
}


#' @title Create a path to a tempory PhyloXML
#' @description An internal function that creates a temporary directory to hold
#' the PhyloXML file for a phylogram. The directory created is only intended to
#' hold the PhyloXML, and no other files, because in order for htmlwidgets to
#' make the PhyloXML available in the HTML it copies the entire parent
#' directory.
#' @return A filepath that can be used to save a PhyloXML.
#' @keywords internal
phyloxml_path <- function() {
  # Create a temp dir for the phyloxml file
  phyloxml_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(phyloxml_tmpdir)  # htmlwidgets copies entire dir to browser
  
  xml_file <- tempfile(pattern="phyloxml-", tmpdir=phyloxml_tmpdir,
                       fileext=".xml")
}


#' @title Convert a Newick file to a PhyloXML
#' @description An internal function that utilizes the Forester Java library to
#' perform the conversion. The converted PhyloXML file will be saved to a
#' tempory directory.
#' @param newick_file A path to the Newick file, likely created by the
#' \code{\link{create_tree}} function.
#' @return A path to the PhyloXML file.
#' @keywords internal
newick_to_phyloxml <- function(newick_file, verbose, verbose_dir) {
  xml_file <- phyloxml_path()
  forester <- system.file("java", "forester_1038.jar", package="receptormarker")
  system(sprintf(
    "java -cp %s org.forester.application.phyloxml_converter -f=nn -ni %s %s",
    forester,
    newick_file,
    xml_file
    ), ignore.stdout=!verbose, ignore.stderr=!verbose
  )
  # Also write the phyloxml to the verbose folder if the user wants it
  if (verbose && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  xml_file
}


#' @title Create a PhyloXML using Biopython
#' @description An internal function that utilizes a custom Python script and
#' Biopython to create a PhyloXML from a multiple sequence alignment.
#' @template -ms_alignment 
#' @param verbose \code{TRUE} or \code{FALSE}, depending on whether or not the
#' stdout and stderr of the Python command should be captured by R. As of now,
#' seems to have no effect.
#' @return A path to the PhyloXML file.
#' @keywords internal
phyloxml_via_biopython <- function(ms_alignment, dist, verbose) {
  xml_file <- phyloxml_path()
  phyloxml_from_msa <- system.file("py", "phyloxml_from_msa.py",
                                   package="receptormarker")
  system(sprintf(
    "python %s --msa %s --dest %s --distance %s",
    phyloxml_from_msa,
    ms_alignment,
    xml_file,
    dist
    ),
    ignore.stdout=!verbose,
    ignore.stderr=!verbose
  )
  xml_file
}


#' @title Remove a dedupe hash from a PhyloXML
#' @description An internal function that looks for a hash in all \code{<name>}
#' nodes of a PhyloXMl and removes the hash and any remaining characters after
#' it. Used in conjuction with the function \code{\link{msa}} when its argument
#' \code{dedupe_hash} is used.
#' @template -xml_file
#' @param hash A string to remove in the PhyloXML. In the XML node(s) in which
#' \code{hash} is found, any characters following the \code{hash} will also be
#' removed.
#' @return A path to the PhyloXML file with the hash removed.
#' @seealso \code{\link{msa}}
#' @keywords internal
remove_phyloxml_hash <- function(xml_file, hash) {
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  named_nodes <- XML::getNodeSet(root, "//ns:name", namespaces=ns)
  # Grab the separating period + hash + rest of chars until end of str
  regex_find <- paste0("\\.", hash, ".*$")  
  lapply(named_nodes, function(n) {
    XML::xmlValue(n) <- gsub(regex_find, "", XML::xmlValue(n))
  })
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}


#' @title Calculate the arc size for a radial phylogram
#' @description An internal function that takes into account the number of
#' sequences in a PhyloXML and adjusts the arc size the achieve uniform spacing
#' between sequences.
#' @template -xml_file
#' @return An integer value representing the arc size for the phylogram
#' @seealso \code{\link{radial_phylo}}
#' @keywords internal
calculate_arc_size <- function(xml_file) {
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
  if (num_elements == 0) {
    stop("Cannot generate phylogram: no sequences to plot", call.=FALSE)
  } else if (num_elements <= 15) {
    arc <- 300
  } else if (num_elements <= 30) {
    arc <- 270
  } else if (num_elements <= 50) {
    arc <- 200
  } else if (num_elements <= 100) {
    arc <- 80
  } else if (num_elements <= 150) {
    arc <- 30
  } else if (num_elements <= 1000) {
    arc <- 20
  } else {
    arc <- 10
  }
  arc
}


#' @title Calculate the canvas size for a radial phylogram
#' @description An internal function that takes into account the length of
#' sequences in a PhyloXML, whether or not the phylogram is condensed, and the
#' number of rings surrounding the phlogram in order to calculate the optimal
#' canvas size for drawing the SVG. The calculations tend to err on the side of
#' creating a larger canvas than a smaller one in order to not cut off parts of
#' the phylogram when it is drawn.
#' @template -xml_file
#' @template -condensed
#' @param rings A user-supplied named character vector representing rings to add
#' to a phylogram.
#' @return An integer value representing both the pixel width and height for the
#' canvas.
#' @seealso \code{\link{radial_phylo}}
#' @keywords internal
calculate_canvas_size <- function(xml_file, condensed, rings) {
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
  # nolint start
  # Similarly: xpathApply(doc, "/ns:phyloxml//ns:name", xmlValue, namespaces=ns)
  # nolint end
  if (num_elements == 0) {
    stop("Cannot generate phylogram: no sequences to plot", call.=FALSE)
  } else if (num_elements <= 30) {
    base_size <- 1450
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
  # Need to add extra room if there are outer rings (~100 pixels per ring works,
  # unless there is a ring that is annotating 'all', then that's going to create
  # many additional rings).
  if (!is.null(rings)) {
    canvas_size <- canvas_size + length(rings) * 120
  }
  if (!(is.null(rings)) && num_elements <= 90) {
    canvas_size <- canvas_size + 200
  }
  if (condensed == TRUE) {
    canvas_size <- canvas_size + 200
  }
  canvas_size
}


#' @title Add bars data to a PhyloXML containing unique sequences
#' @description An internal function that performs the final step of createing a
#' "condensed" phylogram. It calculates sequence frequencies and adds the XML
#' data in order to show bars (representing clonal expansion) and a tooltip on
#' the radial phylogram. 
#' @param xml_file A path to a PhyloXML file that has been built using only
#' unique sequences.
#' @param seqs A character vector containing non-unique sequences. This function
#' will use the frequencies of sequences in this vector to add bars to the
#' corresponding sequences in the PhyloXML. The \code{xml_file} representing the
#' PhyloXML should have been built from the unique sequences in this vector,
#' i.e. with \code{unique(seqs)}.
#' @return A path to the PhyloXML file annotated with bars and tooltip data.
#' @seealso \code{\link{radial_phylo}}
#' @keywords internal
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
  #   <branch_length>0.141994136</branch_length>
  #   <chart>
  #     <content>3</content>  # This makes the bar length
  #   </chart>
  #   <annotation>
  #     <desc>Expansions: 3</desc>  # This is the Unitip tooltip
  #   </annotation>
  # </clade>
  lapply(named_nodes, function(n) {
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
  })
  # Add the bar styles, which need to be appended to the phylogeny
  xpath <- "/ns:phyloxml/ns:phylogeny"  # nolint
  phylogeny_set <- XML::getNodeSet(root, xpath, namespaces=ns)
  phylogeny <- phylogeny_set[[1]]
  # <render>
  #  <charts>
  #    <content fill="#666" type="bar" width="0.6"/>
  #  </charts>
  #  <styles>
  #    <barChart fill="#333" stroke-width="0"/>
  #  </styles>
  # </render>
  render_el <- XML::newXMLNode("render")
  charts_el <- XML::newXMLNode("charts")
  content_el <- XML::newXMLNode("content", attrs=c(fill="#666",
                                                   type="bar",
                                                   width="0.6"))
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


#' @title Index a vector and allow the indexation to wrap if necessary
#' @description An internal function that allows for indexation of a vector to
#' wrap back to the beginning of a vector if the index is greater than the
#' vector's length.
#' @param index An index to use to subset the vector
#' @param v A vector, either named or not named
#' @return A list containing the \code{"name"} and \code{"val"} of the element.
#' @examples
#' \dontrun{
#' color <- c(white="#FFF", black="000", blue="#00F")
#' index_with_wrap(5, color)
#' }
#' @keywords internal
index_with_wrap <- function(index, v) {
  indx <- index %% length(v)
  if (indx == 0 ) indx <- length(v)
  list(name=names(v)[indx], val=v[indx])
}


#' @title Add rings to a PhyloXML for each unique value in a single column
#' @description An internal function that creates an outer ring for each
#' unique value in a single column. Ring annotations are not mutually exclusive,
#' so if a single sequence repeats and it has multiple labels in the annotated
#' column then the single sequence will have multiple rings annotated.
#' @template -xml_file
#' @template -seqs
#' @param d_clean The cleaned data frame from which the \code{seqs} were
#' extracted, both preferably cleaned by \code{\link{clean_data}}.
#' @param seqs_col Either an integer corresponding to a column index or a string
#' corresponding to a column name in d that contains the sequences.
#' @template -rings
#' @param colors A named character vector where the names typically
#' describe a color and the values correspond to 6-character hex color codes
#' (preceded with "#"). This should correspond to the \code{color_wheel}
#' parameter passed to \code{\link{radial_phylo}}.
#' @param ring_thickness An integer representing the thickness of the rings.
#' @template -condensed
#' @return A list containing the path to the PhyloXML \code{xml_file} annotated
#' with rings data, and the \code{ring_map} matching unique values in the
#' annotated column with their colors.
#' @seealso \code{\link{radial_phylo}}, \code{\link{add_phylo_outer_rings}}
#' @keywords internal
add_phylo_outer_rings_all <- function(xml_file, seqs, d_clean, seqs_col,
                                      rings, colors, ring_thickness,
                                      condensed) {
  # Need to create a color map, mapping each unique ring value to a color in
  # the color wheel you're using.
  ring_col <- names(rings)[[1]]
  ring_values <- unique(d_clean[, ring_col])
  ring_values[is.na(ring_values)] <- "N/A (NULL)"
  ring_map <- vapply(seq_along(ring_values), FUN.VALUE=character(1),
                     FUN=function(x) {
    index_with_wrap(x, colors)[["name"]]
  })
  names(ring_map) <- ring_values
  
  # If condensed == TRUE, the sequences in xml_file are unique; in order to get
  # attributes for each ring, for each sequence you have to look at each
  # occurrence of it and annotate appropriately.
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
      if (condensed == FALSE) {
        # Annotate the appropriate ring based on the value in the ring column
        data_row <- seq_data_rows[i, ]
        # Get the value in the current row's ring column
        val <- data_row[[ring_col]]
        if (is.na(val)) {
          val <- "N/A (NULL)"
        }
        # Get the index of that value in the ring map
        idx <- which(names(ring_map) == val)
        # Add the appropriate XML node to color the appropriate ring
        outer_ring_name <- paste0(c("outergroup", idx), collapse="")
        outer_ring_el <- XML::newXMLNode(outer_ring_name, ring_map[[idx]])
        # No <chart> node if condensed == FALSE so add it now
        chart_el <- XML::newXMLNode("chart")
        chart_el <- XML::addChildren(chart_el, kids=c(outer_ring_el))
        parent <- XML::addChildren(parent, chart_el)
      } else if (condensed == TRUE) {
        # For this sequence, get the unique values in the ring column and
        # annotate each of them.
        unique_ring_col_vals <- unique(seq_data_rows[, ring_col])
        unique_ring_col_vals[is.na(unique_ring_col_vals)] <- "N/A (NULL)"
        rings_to_add <- lapply(unique_ring_col_vals, function(val) {
          # Very similar to above
          idx <- which(names(ring_map) == val)
          outer_ring_name <- paste0(c("outergroup", idx), collapse="")
          outer_ring_el <- XML::newXMLNode(outer_ring_name, ring_map[[idx]])
          outer_ring_el
        })
        # The <chart> node exists if condensed == TRUE
        chart_el <- XML::xmlElementsByTagName(parent, "chart")[[1]]
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
  # Add ring styles. If condensed == FALSE there won't be a <styles> or
  # <charts> section yet. Have to create them box. <charts> section dictates
  # how the rings will be drawn and what data type they are. <style> matches
  # up with the colors ('orange', 'blue', etc.)
  xpath <- "/ns:phyloxml/ns:phylogeny"  # nolint
  phylogeny_set <- XML::getNodeSet(root, xpath, namespaces=ns)
  phylogeny <- phylogeny_set[[1]]
  # Outer rings color styles
  outergroups <- vector(mode="list", length=length(ring_map))
  for (q in c(1:length(ring_map))) {
    outergroup_name <- paste0(c("outergroup", q), collapse="")
    outergroups[[q]] <- XML::newXMLNode(outergroup_name, attrs=c(disjointed=q,
                                              type="binary",
                                              thickness=ring_thickness))
  }
  # List names are not preserved in lapply...have to use indices
  color_styles <- lapply(seq_along(ring_map), function(y, n, i) {
    hex_color <- colors[[y[[i]]]]
    XML::newXMLNode(y[[i]], attrs=c(fill=hex_color, opacity="#0.95",
                                    stroke="#DDDDDD"))
  }
  , y=ring_map, n=names(ring_map)
  )
  if (condensed == FALSE) {
    render_el <- XML::newXMLNode("render")
    charts_el <- XML::newXMLNode("charts")
    charts_el <- XML::addChildren(charts_el, kids=outergroups)
    styles_el <- XML::newXMLNode("styles")
    styles_el <- XML::addChildren(styles_el, kids=color_styles)
    render_el <- XML::addChildren(render_el, charts_el, styles_el)
    phylogeny <- XML::addChildren(phylogeny, render_el)
  } else if (condensed == TRUE) {
    xpath <- "/ns:phyloxml/ns:phylogeny/ns:render/ns:charts"  # nolint
    charts_el <- XML::getNodeSet(root, xpath, namespaces=ns)[[1]]
    charts_el <- XML::addChildren(charts_el, kids=outergroups)
    xpath <- "/ns:phyloxml/ns:phylogeny/ns:render/ns:styles"  # nolint
    styles_el <- XML::getNodeSet(root, xpath, namespaces=ns)[[1]]
    styles_el <- XML::addChildren(styles_el, kids=color_styles) 
  }
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
  # Generate better ring color map to return
  color_map <- colors[ring_map]
  names(color_map) <- names(ring_map)
  list("xml_file"=xml_file, "ring_map"=color_map)
}


#' @title Add rings to a PhyloXML corresponding to a matched value
#' @description An internal function that adds a given number of rings to a
#' condensed or not-condensed radial phylogram. Using the sequences in the
#' PhyloXML, it examines the provided data \code{d_clean} and adds the XML data
#' in order to create rings on a radial phylogram whenever the criteria provided
#' by the argument \code{rings} is found to be true in the data. The rings will
#' be colored using a preset color scheme.
#' @template -xml_file
#' @template -seqs
#' @param d_clean The cleaned data frame from which the \code{seqs} were
#' extracted, both preferably cleaned by \code{\link{clean_data}}.
#' @param seqs_col Either an integer corresponding to a column index or a string
#' corresponding to a column name in d that contains the sequences.
#' @template -rings
#' @param colors A named character vector where the names typically
#' describe a color and the values correspond to 6-character hex color codes
#' (preceded with "#"). This should correspond to the \code{color_wheel}
#' parameter passed to \code{\link{radial_phylo}}.
#' @param ring_thickness An integer representing the thickness of the rings.
#' @template -condensed
#' @return A path to the PhyloXML file annotated with rings data.
#' @seealso \code{\link{radial_phylo}}, \code{\link{add_phylo_outer_rings_all}}
#' @keywords internal
add_phylo_outer_rings <- function(xml_file, seqs, d_clean, seqs_col,
                                  rings, colors, ring_thickness, condensed) {
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
  xpath <- "/ns:phyloxml/ns:phylogeny"  # nolint
  phylogeny_set <- XML::getNodeSet(root, xpath, namespaces=ns)
  phylogeny <- phylogeny_set[[1]]
  # Outer rings color styles
  outergroups <- vector(mode="list", length=length(rings))
  for (q in c(1:length(rings))) {
    outergroup_name <- paste0(c("outergroup", q), collapse="")
    outergroups[[q]] <- XML::newXMLNode(outergroup_name, attrs=c(disjointed=q,
                                                     type="binary",
                                                     thickness=ring_thickness))
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
    xpath <- "/ns:phyloxml/ns:phylogeny/ns:render/ns:charts"  # nolint
    charts_el <- XML::getNodeSet(root, xpath, namespaces=ns)[[1]]
    charts_el <- XML::addChildren(charts_el, kids=outergroups)
    xpath <- "/ns:phyloxml/ns:phylogeny/ns:render/ns:styles"  # nolint
    styles_el <- XML::getNodeSet(root, xpath, namespaces=ns)[[1]]
    styles_el <- XML::addChildren(styles_el, kids=color_styles)
  }
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}

#' @title Remove sequence labels from PhyloXML
#' @description An internal function that removes the labels from the phylogram.
#' @template -xml_file
#' @return A path to the PhyloXML file with sequences removed.
#' @seealso \code{\link{radial_phylo}}
#' @keywords internal
remove_phyloxml_labels <- function(xml_file) {
  doc <- XML::xmlParse(xml_file)
  root <- XML::xmlRoot(doc)
  ns <- c(ns=XML::xmlNamespace(root))
  # Get all elements that don't start with Biopython's "InnerXXX"
  named_nodes <- XML::getNodeSet(root,
                                 "//ns:name[not(starts-with(text(), 'Inner'))]",
                                 namespaces=ns)
  
  # Remove all sequences from <name> nodes
  lapply(named_nodes, function(n) {
    XML::xmlValue(n) <- ""
  })
  
  XML::saveXML(doc=doc, file=xml_file,
               prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
               indent=FALSE)
  xml_file
}

#' @title Create an interactive radial phylogram
#' @description Creates a JavaScript-based radial phylogram in RStudio or a
#' browser window. Intended for use with relatively short amino acid sequences,
#' particularly those representing the variable chains or complementary
#' determining regions (CDRs) of antibody or T-cell receptors. Provides several
#' parameters for customizing the phylogram by visually representing clonal
#' expansion and/or visually representing phenotypic traits of individual cells
#' or sequence populations. The phylogram is made available for download as a
#' PNG.
#' @details It is suggested to have Biopython installed and available to R when
#' using this function. If Biopython is not available, a warning will be output
#' when loading the package and running this function. Biopython is able to
#' create a more visually appealing and "smarter" tree using the UPGMA
#' algorithm. If Biopython is not available, \code{\link[ape]{bionj}} will be
#' used to create a tree purely in R.
#' 
#' It is important to note the difference between \code{canvas_size} and
#' \code{scale} when adjusting those arguments. The phylogram is drawn on an SVG
#' canvas. As more sequences, rings, or bars are added to the phylogram, the
#' radius of the phylogram must increase in order to prevent all that
#' information from overlapping. Consequently, the SVG canvas needs to grow to
#' accommodate the growing phylogram. If the growth of the canvas does not keep
#' up with the growth of the phylogram, the phylogram will either be cut off at
#' the horizontal and vertical edges of the canvas (which is akin to an HTML
#' \code{<div>}) or squashed by the small canvas. The argument
#' \code{canvas_size} has everything to do with the initial draw of the
#' phylogram. In contrast, \code{scale} has nothing to do with the size of the
#' phylogram during its initial draw, but rather everything to do with the
#' scaling of the phylogram \emph{after it has already been drawn}. The
#' phylogram will only be drawn once, since scaling does not initiate a redraw.
#' Therefore, any anomolies in the drawing of the phylogram should be addressed
#' using \code{canvas_size}.
#'
#' \strong{TIP:} When experimenting or crafting a large or complex phylogram,
#' such as a condensed phylogram with several annotations, set \code{fast=TRUE}
#' to iterate quickly and when the ideal settings have been found change to
#' \code{fast=FALSE} to build the final phylogram.
#'
#' \strong{TIP:} If the desired outcome is a PNG image of the phylogram, set
#' \code{scale=FALSE, browser=TRUE} and then use the "Download PNG" button in
#' the upper-left of the browser window to save the phylogram as a full-sized
#' PNG.
#' @template -d
#' @template -seqs_col
#' @param dist A string representing a distance matrix to use for calculating
#' distances from the sequence alignment. Only applicable if Biopython is
#' installed and available. Must be one of "benner6", "benner22",
#" "benner74", "blosum100", "blosum30", "blosum35", "blosum40", "blosum45",
#" "blosum50", "blosum55", "blosum60", "blosum62", "blosum65", "blosum70",
#" "blosum75", "blosum80", "blosum85", "blosum90", "blosum95", "feng", "fitch",
#" "genetic", "gonnet", "grant", "ident", "johnson", "levin", "mclach",
#" "miyata", "nwsgappep", "pam120", "pam180", "pam250", "pam30", "pam300",
#" "pam60", "pam90", "rao", "risler", "structure".
#' @param condense \code{TRUE} or \code{FALSE}, depending on whether or not the
#' radial phylogram should only contain unique sequences. If \code{TRUE}, clonal
#' expansion (i.e. sequence frequency data) will be represented by orthogonal
#' vertical bars on the perimeter of the phylogram.
#' @template -rings
#' @param canvas_size An integer representing the desired width and height of
#' the SVG canvas in pixels. Defaults to \code{"auto"} to automatically
#' estimate the appropriate SVG canvas size. It is suggested to leave this at
#' the default value and then adjust only if necessary. If parts of the
#' phylogram are cut off, \code{canvas_size} should be increased. See the
#' Details section for more information.
#' @param font_size An integer font size for the phylogram's labels. It is
#' suggested to leave this at the default value and then adjust only if
#' necessary.
#' @param tree_margin A number between 0 and 1 indicating the size of the
#' tree with respect to the phylogram. To give more focus on the
#' sequence labels, increase this value and increase \code{font_size}.
#' @param arc An integer betwen 1 and 359 indicating the size of the split in
#' the phylogram, in degrees. Defaults to \code{"auto"} to automatically
#' estimate the appropriate arc size.
#' @param color_wheel A named character vector where the names typically
#' describe a color and the values correspond to 6-character hex color codes
#' (preceded with "#"). These colors are used for the outer \code{rings} of the
#' phylogram. If more colors than those supplied are needed, the colors in the
#' wheel will be recycled, starting from the beginning.
#' @param ring_thickness An integer representing the thickness of the rings.
#' @param scale \code{TRUE} or \code{FALSE}, depending on whether or not the
#' phylogram should scale to fit the browser or RStudio Viewer window. If
#' \code{FALSE} and the phylogram is on a large canvas, it will be necessary to
#' scroll to see the entire canvas.
#' @param label \code{TRUE} or \code{FALSE}, depending on whether or not
#' sequence labels should be present in the phylogram.
#' @template -browser
#' @param verbose \code{TRUE} or \code{FALSE}. If \code{TRUE} additional output
#' is printed to the R console and the sequences, multiple sequence alignment,
#' and PhyloXML are written to a folder in the working directory.
#' @param fast \code{TRUE} or \code{FALSE}, depending on whether or not the fast
#' NJ (neighbor joining) algorithm in R should be used to build the phylogram.
#' If \code{FALSE} and Biopython is installed, Biopython's UPGMA algorithm will
#' be used instead. This method is slower, but it will generate "smarter" trees
#' and therefore is the default setting. Only applicable if Biopython is
#' installed and available.
#' @import htmlwidgets
#' @examples
#' # Create a condensed radial phylogram with outer-ring annotations
#' data(tcr)  # Packaged data set, a data.frame from a CSV file
#' radial_phylo(tcr, 'seqs', condense=TRUE, rings=c(FOXP3=1, GATA3=1),
#' browser=TRUE)
#' @export
radial_phylo <- function(d, seqs_col=NULL, dist="ident", condense=FALSE,
                         rings=NULL, canvas_size="auto", font_size=12,
                         tree_margin=0.28, arc="auto",
                         color_wheel=c(green="#82A538", orange="#B1903C",
                                       blue="#626DBC", pink="#B5598F",
                                       turqoise="#0DA765", yellow="#E9E700",
                                       purple="#64059C", tangerine="#DE531C"),
                         ring_thickness=10,
                         scale=TRUE, label=TRUE, browser=FALSE, verbose=FALSE,
                         fast=FALSE) {
  
  check_muscle(level="stop")
  biopy_existence <- check_bio_python(level="warn")
  
  # Validate function parameters
  validate_not_null(list(d=d, dist=dist, condense=condense,
                         canvas_size=canvas_size, font_size=font_size,
                         tree_margin=tree_margin, arc=arc,
                         color_wheel=color_wheel, ring_thickness,
                         scale=scale, label=label, browser=browser,
                         verbose=verbose, fast=fast))
  validate_canvas_size(canvas_size)
  validate_font_size(font_size)
  validate_tree_margin(tree_margin)
  validate_ring_thickness(ring_thickness)
  validate_arc(arc)
  validate_true_false(list(condense=condense, scale=scale, browser=browser,
                           label=label, verbose=verbose, fast=fast))
  validate_distance(dist)
  validate_d_seqs(d, seqs_col)
  validate_rings(rings, d)
  validate_color_wheel(color_wheel)
  
  ring_map <- ""  # Used to make color legend for ring annotations if they exist
  
  # Not necessary to have as func parameters; these will get set automatically
  width <- NULL
  height <- NULL
  
  # Create verbose dir
  if (verbose) {
    verbose_dir <- tempfile("radial_phylo-", tmpdir=getwd(), fileext="")
    dir.create(verbose_dir)
  }
  
  # Step 1: Clean the user-supplied data and the sequences
  clean <- clean_data(d, seqs_col, verbose, verbose_dir, verbose_format="fasta")
  seqs <- clean[["seqs"]]
  if (condense == TRUE) {
    seqs <- unique(seqs)
  }
  # Step 2: Do a multiple sequence alignment (MSA)
  dedupe_hash <- paste0(sample(c(0:9, letters), 10, replace=TRUE), collapse="")
  ms_alignment <- msa(seqs, verbose, verbose_dir, dedupe_hash)
  if (is.null(biopy_existence) || fast == TRUE) {
    # Step 3: Compute pairwise distances of the aligned sequences
    dist_matrix <- compute_dist_matrix(ms_alignment[["as_string"]], seqs)
    # Step 4: Calculate a distance tree and write it as .newick
    newick_file <- create_tree(dist_matrix)
    # Step 5: Convert the .newick to phylo.xml
    xml_file <- newick_to_phyloxml(newick_file, verbose, verbose_dir)
  } else {
    xml_file_with_hash <- phyloxml_via_biopython(ms_alignment[["file"]], dist,
                                                 verbose)
    xml_file <- remove_phyloxml_hash(xml_file_with_hash, dedupe_hash)
  }
  
  # Condense phylo if condense == TRUE
  if (condense == TRUE) {
    # Send the non-unique sequences to be able to count how many there are
    xml_file <- add_bars_to_condensed_phyloxml(xml_file, clean[["seqs"]])
  }
  
  # Add outer-rings to phylo if there are any requested
  if (!is.null(rings) && "all" %in% rings) {
    r <- add_phylo_outer_rings_all(xml_file, seqs, clean[["d"]],
                                   seqs_col, rings, color_wheel, ring_thickness,
                                   condense)
    xml_file <- r[["xml_file"]]
    ring_map <- r[["ring_map"]]
  } else if (!is.null(rings)) {
    xml_file <- add_phylo_outer_rings(xml_file, seqs, clean[["d"]], seqs_col,
                                      rings, color_wheel, ring_thickness,
                                      condense)
  }
  
  # Calculate canvas size based on number of nodes in phylo.xml. 
  if (canvas_size == "auto") {
    # Normally canvas size is calculated by the length of the 'rings' variable,
    # but when you allow all values of a column to be annotated with col='all',
    # the number of rings get expanded into the number of unique values in that
    # column. So, need to account for that.
    total_rings <- c()
    for (ring_name in names(rings)) {
      ring_val <- rings[[ring_name]]
      if (ring_val == "all") {
        expanded_rings <- unique(clean[["d"]][, ring_name])
        total_rings <- append(total_rings, expanded_rings)
      } else {
        total_rings <- append(total_rings, c(ring_name=ring_val))
      }
    }
    canvas_size <- calculate_canvas_size(xml_file, condense, total_rings) 
  }
  
  if (arc == "auto") {
    arc <- calculate_arc_size(xml_file)
  }
  
  if (!label) {
    xml_file <- remove_phyloxml_labels(xml_file) 
  }
  
  # Also write the phyloxml to the verbose folder if the user wants it
  if (verbose && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir)
  }
  
  # Forward options to radial_phylo.js using 'x'
  x <- list(
    canvas_size = canvas_size,
    scale = scale,
    font_size = font_size,
    tree_margin = tree_margin,
    arc = arc,
    legend_values = names(ring_map),
    legend_colors = unname(ring_map)
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


#' @title Shiny bindings for radial phylograms
#' @description Output and render functions for using \code{\link{radial_phylo}}
#' within Shiny applications and interactive Rmd documents.
#' @param outputId The output variable to read from.
#' @param width,height A valid CSS unit, such as \code{"100\%"},
#' \code{"600px"}, \code{"auto"}, or a number (which will be coerced to a string
#' and have \code{"px"} appended to it).
#' @param expr An expression that generates a radial phylogram.
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#' is useful if you want to save an expression in a variable.
#' @name radial_phylo_shiny
#' @export
# nolint start
radial_phyloOutput <- function(outputId, width = "100%", height = "400px"){
  shinyWidgetOutput(outputId, "radial_phylo", width, height,
                    package="receptormarker")
# nolint end
}


#' @rdname radial_phylo_shiny
#' @export
# nolint start
renderRadial_phylo <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr) # force quoted
  shinyRenderWidget(expr, radial_phyloOutput, env, quoted = TRUE)
# nolint end
}
