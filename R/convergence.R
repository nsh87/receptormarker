# Minglu Ma
# This script compare gene and find the similary parts
# Date: 11/08/2015


# prepare for convergence
# Minglu
# For test
string <- "This is my string string"
d <- strsplit(string, " ")[[1]]

# clean the input vector to unique string and save to txt file
clean_input <- function(d) {
      print(d)
      d_unique <- unique(d)
      print(d_unique)
      write(d_unique, "/Users/mingluma/2015Fall/receptormarker/output.txt", sep="\n")
}

# test system function  
x=system('perl -e “print 2 + 4″', intern=TRUE)

# run perl script
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



  
