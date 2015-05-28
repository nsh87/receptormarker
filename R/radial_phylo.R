#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
# allow users to set viewer.suppress to FALSE to see the thing in RStudio
radial_phylo <- function(dataFile, radius='fill', autoResize=FALSE,
                         suppressViewer=FALSE, width=NULL, height=NULL) {

  u <- read.csv(dataFile)
  seqscol <- 3
  seqs <- u[u[, seqscol] != '', ]
  w <- with(seqs, paste0('>', seqs[, seqscol], '\n', seqs[, seqscol], collapse='\n'))
  seqsFile <- tempfile(pattern='seqs-', tmpdir=getwd(), fileext='.txt')
  write(w, seqsFile)
  
  alignedFile <- tempfile(pattern='aligned-', tmpdir=getwd(), fileext='.txt')
  s <- as.character(seqs[, seqscol])
  ss <- Biostrings::AAStringSet(s)
  alignment <- muscle::muscle(stringset=ss, quiet=TRUE)
  caretseqs <- as.character(lapply(s, function(x) { paste0('>', x) }))
  fast <- lapply(seq_along(caretseqs), function(y, i) {
    caret <- y[i]
    algn <- as.character(alignment)[i]
    c(caret, algn)
  }, y=caretseqs)
  faster <- unlist(fast)
  write(faster, alignedFile)
  
  aligned <- seqinr::read.alignment(alignedFile, format='fasta')
  aligned_dist <- as.matrix(seqinr::dist.alignment(x=aligned, matrix='identity'))
  
  tree <- ape::nj(aligned_dist)
  new_tree <- ape::as.phylo(tree)
  newickFile <- tempfile(pattern='tree-', tmpdir=getwd(), fileext='.newick')
  ape::write.tree(phy=new_tree, file=newickFile)
  
  xmldir = paste0(getwd(), '/phyloxml')
  dir.create(xmldir)
  xmlFile <- tempfile(pattern='phyloxml-', tmpdir=xmldir, fileext='.xml')
  forester <- system.file("java", "forester_1038.jar", package="receptormarker")
  system(sprintf("java -cp %s org.forester.application.phyloxml_converter -f=nn -ni %s %s", forester, newickFile, xmlFile))
  
  # Determine what the parameter 'radius' is
  err = "The argument 'radius' is invalid"
  radius_options <- c('fill')
  tryCatch({
    if (!is.element(radius, radius_options) && (radius != floor(radius))) {
      stop(err, call.=FALSE)
    }
  },
  error = function(e) {
    stop(err, call.=FALSE)
  }
  )
  
  # Warn user if setting 'radius' and autoResize=TRUE at same time
  warn = "Using autoResize=TRUE will override the 'radius' parameter"
  tryCatch({
    if (is.numeric(radius) && autoResize == TRUE) {
      warning(warn, call.=FALSE)
    }
  })
  
  # forward options using x
  x = list(
    radius = radius,
    autoResize = autoResize
  )
  
#   xmlFileSplit <- strsplit(xmlFile, '/', fixed=TRUE)
#   xmlFileNameOnly <- xmlFileSplit[[1]][length(xmlFileSplit[[1]])]
#   
#   
   phyloxml <- htmltools::htmlDependency(
     name = 'phyloxml',
     version = '1.0',
#     src = system.file('htmlwidgets/lib/phyloxmls', package='receptormarker'),
     src = c(file=dirname(xmlFile)),
     attachment = list(xml=basename(xmlFile))
   )
#   dep <- htmltools::copyDependencyToDir(phyloxml, outputDir='htmlwidgets/lib')

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
