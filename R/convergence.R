# For test save the input vector to unique string and save to txt file
#input_d <- read.csv("/Users/mingluma/2015Fall/receptormarker/vdjfasta/bin/sample.txt", sep = " ", header = FALSE)

# return the temp directory
save_input <- function(d) {
  # print(d)
  d_unique <- unique(d)
  # print(d_unique)
  input_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  input <- tempfile(pattern="sample_", tmpdir=input_tmpdir, fileext=".txt")
  write(d_unique, input, sep="\n")
  input_tmpdir
}

# run perl script isd, TODO:
run_perl <- function(input_file){
  #set the path of hmm and blast to system PATH
  #written by Daniel Zhang till "##"
  sys_path <- Sys.getenv("PATH")
  hmm_tool <- system.file("tools/hmmer-3.1b2/binaries", package = "receptormarker")
  blstn_tool <- system.file("tools/ncbi/blast/bin", package = "receptormarker")
  Sys.setenv( PATH = paste(sys_path,hmm_tool,blstn_tool,sep=":"))
  ##
  perl_file <- system.file("perl/vdjfasta/bin", "convergence-pipeline.pl",
                           package="receptormarker")
  system(sprintf("perl %s --textfile=%s",
                 perl_file,
                 input_file
  )
  # verbose means?
  # ignore.stdout=!verbose, 
  # ignore.stderr=!verbose
  )
}

# get xml saved path
convergencexml_path <- function() {
  # Create a temp dir for the convergencexml file
  convergencexml_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  dir.create(convergencexml_tmpdir)  # htmlwidgets copies entire dir to browser 
  xml_file <- tempfile(pattern="convergencexml-", tmpdir=convergencexml_tmpdir,
                       fileext=".xml")
}

#read output and write to xml file as the requirement from siwei.
read_output <- function(path) {
  xml_file <- convergencexml_path()
  test <- read.csv(path, sep = '\t', header = FALSE)
  
  # For test, I hard code the output path
  # test <- read.csv("/Users/mingluma/2015Fall/receptormarker/vdjfasta/bin/sample-convergence-groups.txt", 
  # sep = '\t', header = FALSE)
  test <- as.data.frame(test)
  # Get the row with longest nodes, just read the 11st row
  ############################# ONLY READ ONE ROW ##################################
  # TODO： Should be modify to loop all the rows
  nodes_list <- test[14,][3]
  node_char <- as.character(nodes_list$V3)
  node_char <- strsplit(node_char, " ")
  df_nodes <- as.data.frame(node_char)
  #dim(df_nodes) each seq is df_node[i,], i from 1 to rownums
  # if I get df_nodes as input. I will have
  rownums<-nrow(df_nodes)
  
  # Generate XML Tree
  node = XML::newXMLNode("graphml")
  keynode = XML::newXMLNode("key", 
                       attrs = c(id="label", "for"="all", "attr.name"="label", 
                                 "attr.type"="string"), parent = node)
  graphnode = XML::newXMLNode("graph", attrs = c(id="0", 
                                            "edgedefault"="undirected"), parent = node)
  kidsnode = lapply(c(1:rownums),
                    function(x)
                      XML::newXMLNode("node", attrs = c(id = x), 
                                 .children = sapply(x, function(x) 
                                   XML::newXMLNode("data", df_nodes[x,],
                                              attrs = c(key ="label"))) ))
  
  XML::addChildren(graphnode, kidsnode)
  if (rownums>1){
    for(i in 1:(rownums-1)){
      kidsedge = lapply(c((i+1):rownums),
                        function(x)
                          XML::newXMLNode("edge","", 
                                     attrs = c("source"= i , "target"=x)))
      XML::addChildren(graphnode, kidsedge)
    }
  }
  #cat(saveXML(node)) ? TODO:check
  XML::saveXML(node, file=xml_file,
          prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
          indent=FALSE)
  xml_file
  file.copy(xml_file, getwd())
}


#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#' @example 
#' 
#' @export
convergence <- function(d, seqs_col=NULL, condense=FALSE, rings=NULL,
                        canvas_size="auto", font_size=12, scale=TRUE,
                        browser=FALSE, verbose=FALSE, fast=FALSE) {
  validate_not_null(list(d=d, condense=condense, canvas_size=canvas_size,
                         font_size=font_size, scale=scale, browser=browser,
                         verbose=verbose, fast=fast))  
  validate_true_false(list(condense=condense, scale=scale, browser=browser,
                           verbose=verbose, fast=fast))  
  validate_sequences(d)
  
  # clean data to unique sequence
  d_unique <- unique(d)
  
  # create temp dir for input and out put
  convergence_tmpdir <- tempfile("", tmpdir=tempdir(), fileext="")
  
  dir.create(convergence_tmpdir)
  
  input <- tempfile(pattern="sample", tmpdir=convergence_tmpdir, fileext=".txt")
  
  #write(d_unique, input, sep="\n")
  write(as.character(d_unique[,1]), input, sep="\n")
  
  
  #get output path
  output <- gsub(pattern = ".txt$", replacement = "-convergence-groups.txt", x = input, ignore.case = T)
  #print(output)
  
  # run perl script
  print("before run_perl")
  run_perl(input)
  print("after run perl")
  
  # Read output file and save to xml
  xml_file <-read_output(output)
  #print(xml_file)
  #################TODO: follow code need to modify? #################
  # Forward options to radial_phylo.js using 'x'
  x <- list(
    canvas_size = canvas_size,
    scale = scale,
    font_size = font_size
  )
  
  
  # Add the phyloxml as an HTML dependency so it can get loaded in the browser
  convergencexml <- htmltools::htmlDependency(
    name = "convergencexml",
    version = "1.0",
    src = c(file=dirname(xml_file)),
    attachment = list(xml=basename(xml_file))
  )
  width=NULL
  height=NULL
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
    dependencies = convergencexml
  )
}

#' Widget output function for use in Shiny
#'
#' @export
convergenceOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'convergence', width, height, package = 'receptormarker')
}

#' Widget render function for use in Shiny
#'
#' @export
renderConvergence <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, convergenceOutput, env, quoted = TRUE)
}
