#!/usr/bin/env sh
log_file="travis/travis_build.log"

# Run custom commands, such as running your test cases
eval "Rscript ./travis/travis_build.R | tee -a '$log_file'"

# Search for errors in the output of your custom commands
err1="ERROR"
err2="WARNING"
err3="Failure (at"
err4="Failure(@"
err5="Error: "
if ! grep -q "$err1\|$err2\|$err3\|$err4\|$err5" $log_file; then
    echo "No errors, warnings, or failures found."
else 
    printf "\n"
    echo "*** grep results **********************"
    grep -n "$err1" $log_file
    grep -n "$err2" $log_file
    grep -n "$err3" $log_file
    grep -n "$err4" $log_file
    grep -n "$err5" $log_file
    echo "ERROR, WARNING, or Failure found. See grep results above."
    printf "\n"
    exit 1
fi
