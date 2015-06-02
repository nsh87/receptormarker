#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(df, seqsCol, radius='fill', fontSize='auto',
                         autoResize=FALSE, suppressViewer=FALSE, width=NULL,
                         height=NULL, verbose=FALSE) {

  # Determine what the parameter 'radius' is
  err = "The argument 'radius' is invalid"
  radius_options <- c('fill')
  tryCatch({
    if (!is.element(radius, radius_options) && (radius != floor(radius)) && length(radius) != 1) {
      stop(err, call.=FALSE)
    } else if (radius < 1) {
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(err, call.=FALSE)
  }
  )
  
  # Get the sequences column from the data.frame
  seqsColErr <- "The argument 'df' and/or 'seqsCol' is invalid" 
  tryCatch({
    if (typeof(seqsCol) == 'character' && length(names(df)) > 0 && length(seqsCol) == 1) {
      seqs <- as.character(df[, seqsCol])
    } else if (seqsCol == floor(seqsCol) && length(seqsCol) == 1) {
      seqs <- as.character(df[, seqsCol])
    } else {
      stop(seqsColErr, call.=FALSE)
    }
  },
  error = function(e) {
    stop(seqsColErr, call.=FALSE)
  }
  )
  
  # Make sure sequences are only alpha characters
  seqsColErr <- "Sequences must only contain characters from A-Z"
  g <- grepl("[^A-Za-z]", as.character(seqs))
  if (sum(g) > 0) {
    stop(seqsColErr, cal.=FALSE)
  }
  
  # Create temp dirs to hold intermediate files unless user want to see them
  tmp_dir <- tempdir()  # Different for each new R session
  # PhyloXML's should go in their own sub-dir because htmlwidgets will copy its 
  # entire parent dir to make the file available in the browser
  phyloxml_tmpdir <- tempfile('', tmpdir=tempdir(), fileext='')
  dir.create(phyloxml_tmpdir)
  # Create verbose dir if user wants the intermediate files
  if (verbose == TRUE) {
    verbose_dir <- tempfile('radial_phylo-', tmpdir=getwd(), fileext='')
    dir.create(verbose_dir)
  }
  
  # Step 1: Clean the data.frame and get the cleaned sequences
  df_clean <- df[df[, seqsCol] != '', ]  # Remove rows with no sequences
  df_clean <- df_clean[complete.cases(df_clean[, seqsCol]), ]  # Remove rows with NA seqs
  seqs <- as.character(df_clean[, seqsCol])
  # Write the sequences if the user wants them
  if (verbose == TRUE) {
    seqs_fasta <- with(df_clean,
                       paste0(">", df_clean[, seqsCol], "\n",
                              df_clean[, seqsCol],
                              collapse="\n"
                              )
    )
    seqs_file <- tempfile(pattern='sequences-', tmpdir=verbose_dir, fileext='.txt')
    write(seqs_fasta, seqs_file)
  }
  
  # Step 2: Do a multiple sequence alignment (MSA)
  seqs_biostring <- Biostrings::AAStringSet(seqs)
  names(seqs_biostring) <- seqs
  if (verbose == TRUE) { print("MUSCLE multiple sequence alignment:") }
  ms_alignment <- muscle::muscle(stringset=seqs_biostring, quiet=!verbose)
  # Write the alignment if the user wants it
  if (verbose == TRUE) {
    aa_str_set <- as(ms_alignment, "AAStringSet")
    msa_file <- tempfile(pattern='msa-', tmpdir=verbose_dir, fileext='.fasta')
    Biostrings::writeXStringSet(aa_str_set, file=msa_file)
  }
  
  # Step 3: Compute pairwise distances of the aligned sequences
  # Need to turn the MSA into class 'alignment'
  alignment <- seqinr::as.alignment(nb=nrow(ms_alignment),
                            nam=seqs,
                            seq=as.character(ms_alignment, use.names=FALSE)
  ) 
  dist_matrix <- as.matrix(seqinr::dist.alignment(x=alignment, matrix='identity'))
  
  # Step 4: Calculate a distance tree and write it as .newick
  dist_tree <- ape::bionj(dist_matrix)
  phylo_tree <- ape::as.phylo(dist_tree)
  newick_file <- tempfile(pattern='tree-', tmpdir=tmp_dir, fileext='.newick')
  ape::write.tree(phy=phylo_tree, file=newick_file)
  # Copy the newick file from the tmp dir to the verbose dir if the user wants it
  if (verbose == TRUE && file.exists(newick_file)) {
    file.copy(newick_file, verbose_dir, copy.date=TRUE)
  }
  
  # Step 5: Convert the .newick to phylo.xml
  xml_file <- tempfile(pattern='phyloxml-', tmpdir=phyloxml_tmpdir, fileext='.xml')
  forester <- system.file("java", "forester_1038.jar", package="receptormarker")
  system(sprintf("java -cp %s org.forester.application.phyloxml_converter -f=nn -ni %s %s", forester, newick_file, xml_file), ignore.stdout=!verbose, ignore.stderr=!verbose)
  # Also write it to the verbose folder if the user wants it
  if (verbose == TRUE && file.exists(xml_file)) {
    file.copy(xml_file, verbose_dir, copy.date=TRUE)
  }
  
  # forward options to radial_phylo.js using 'x'
  x = list(
    radius = radius,
    autoResize = autoResize
  )
  
  # add the phyloxml as an HTML dependency so it can get loaded in the browser
   phyloxml <- htmltools::htmlDependency(
     name = 'phyloxml',
     version = '1.0',
     src = c(file=dirname(xml_file)),
     attachment = list(xml=basename(xml_file))
   )

  # create widget
  htmlwidgets::createWidget(
    name = 'radial_phylo',
    x,
    width = width,
    height = height,
    htmlwidgets::sizingPolicy(
      viewer.padding = 0,
      viewer.suppress = suppressViewer,
      browser.fill = TRUE
    ),
    package = 'receptormarker',
    dependencies = phyloxml
  )
  
}

#' Widget output function for use in Shiny
#'
#' @export
radial_phyloOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'radial_phylo', width, height, package = 'receptormarker')
}

#' Widget render function for use in Shiny
#'
#' @export
renderRadial_phylo <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, radial_phyloOutput, env, quoted = TRUE)
}
