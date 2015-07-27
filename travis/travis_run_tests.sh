#!/usr/bin/env sh
log_file="travis/travis_build.log"

echo "Running package's unit tests"

# Run custom commands, such as running your test cases
eval "Rscript ./travis/travis_run_tests.R | tee -a '$log_file'"

# Search for errors in the output of your custom commands
err1="ERROR"
err2="WARNING"
err3="Failure (at"
err4="Failure(@"
err5="Error: "
err6="Execution halted"
if ! grep -q "$err1\|$err2\|$err3\|$err4\|$err5\|$err6" $log_file; then
    echo "No errors, warnings, or failures found when running unit tests."
else 
    printf "\n"
    echo "*** grep results from unit tests ********************"
    grep -n "$err1" $log_file
    grep -n "$err2" $log_file
    grep -n "$err3" $log_file
    grep -n "$err4" $log_file
    grep -n "$err5" $log_file
    grep -n "$err6" $log_file
    echo "ERROR, WARNING, or Failure found in unit tests. See grep results above."
    printf "\n"
    exit 1
fi
