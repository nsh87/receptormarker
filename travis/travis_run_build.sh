#!/usr/bin/env sh
log_file="travis/travis_build.log"
#R CMD BATCH ./travis_build.R $log_file 
#cat $log_file

# touch $log_file
# 
# # run tail -f in background
# tail -f $log_file # > out 2>&1 &
# 
# # process id of tail command
# tailpid=$!
# 
# # run build
# R CMD BATCH ./travis_build.R $log_file 
# 
# # now kill the tail process
# kill $tailpid

# run checks
eval "Rscript ./travis/travis_build.R | tee -a '$log_file'"

# search for errors
err1="ERROR"
err2="WARNING"
err3="Failure (at"
err4="Failure(@"
if ! grep -q "$err1\|$err2\|$err3\|$err4" $log_file; then
    echo "No errors, warnings, or failures found."
else 
    echo "*** grep results **********************"
    grep -n "$err1" $log_file
    grep -n "$err2" $log_file
    grep -n "$err3" $log_file
    grep -n "$err4" $log_file
    echo "ERROR, WARNING, or Failure found. See grep results above."
    exit 1
fi
