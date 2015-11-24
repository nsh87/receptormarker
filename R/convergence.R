# Minglu Ma
# This script compare gene and find the similary parts
# Date: 11/08/2015

# prepare for convergence
# Minglu
# For test
#string <- "This is my string string"
#d <- strsplit(string, " ")[[1]]

## Copy from radial_phylo.r
validate_sequences <- function(seqs) {
  # Make sure sequences are only alpha characters
  seqs_col_err <- "Sequences must only contain characters from A-Z and a-z"
  g <- grepl("[^A-Za-z]", as.character(seqs))
  if (sum(g) > 0) {
    stop(seqs_col_err, call.=FALSE)
  }
}

# save the input vector to unique string and save to txt file
# return the temp directory
save_input <- function(d) {
      # print(d)
      d_unique <- unique(d)
      # print(d_unique)
      input_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
      input <- tempfile(pattern="sample_", tmpdir=input_tmpdir, fileext=".txt")
      write(d_unique, input, sep="\n")
      return input_tmpdir
}

# test system function  
# x=system('perl -e “print 2 + 4″', intern=TRUE)

# run perl script isd
run_perl <- function(input_file){
  #TODO: how to get output path?
  perl_file <- system.file("pl", "convergence-pipeline.pl",
                           package="receptormarker")
  system(sprintf("perl %s --textfil=%s",
                 perl_file,
                 input_file
                 ),
         ignore.stdout=!verbose, 
         ignore.stderr=!verbose
  )
}

#TODO: After read output, reformat and send to siwei
read_output <- function(path) {
  mydata <- read.table(path)
  print(mydata)
}


# The main function, TODO: understand and clear the function parameters
convergence <- function(d, seqs_col=NULL, condense=FALSE, rings=NULL,
                         canvas_size="auto", font_size=12, scale=TRUE,
                         browser=FALSE, verbose=FALSE, fast=FALSE) {
  validate_not_null(list(d=d, condense=condense, canvas_size=canvas_size,
                         font_size=font_size, scale=scale, browser=browser,
                         verbose=verbose, fast=fast))  
  validate_true_false(list(condense=condense, scale=scale, browser=browser,
                        verbose=verbose, fast=fast))  
  validate_d_seqs(d, seqs_col)

  # clean data to unique sequence
  d_unique <- unique(d)

  # create temp dir for input and out put
  convergence_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(convergence_tmpdir)
  input <- tempfile(pattern="sample", tmpdir=convergence_tmpdir, fileext=".txt")
  write(d_unique, input, sep="\n")

  #?????? TODO slipt()
  output <- tempfile(pattern="sample-convergence-groups", tmpdir=convergence_tmpdir, fileext=".txt")
  # or TODO
  output <- file.path(, convergence_tmpdir)

  # run perl script
  run_perl(input)




#################TODO: follow code need to modify #################
  # Forward options to radial_phylo.js using 'x'
  x <- list(
    canvas_size = canvas_size,
    scale = scale,
    font_size = font_size
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


# devtools::load_all()
  
