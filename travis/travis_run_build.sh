#!/usr/bin/env sh

# Build and check the package (not --as-cran, since they don't allow binaries
# to be included and we include BLAST and HMMER binares)
echo "Building the package with: R CMD build"
R CMD build .
export PKG_TARBALL=$(Rscript -e 'pkg <- devtools::as.package("."); cat(paste0(pkg$package, "_", pkg$version, ".tar.gz"));')

# Check the package
echo "Checking the package with: R CMD check ${PKG_TARBALL}"
R CMD check "${PKG_TARBALL}"

# search for errors in the output of the package check log
export RCHECK_DIR=$(Rscript -e 'cat(paste0(devtools::as.package(".")$package, ".Rcheck"))')
err1="ERROR"
err2="WARNING"
err3="Failure (at"
err4="Failure(@"
err5="Error: "
if ! grep -q -R "$err1\|$err2\|$err3\|$err4\|$err5" "${RCHECK_DIR}/00check.log"; then
    echo "No errors or warnings found when checking the package"
else
    printf "\n"
    echo "*** grep results from package check ********************"
    grep -n "$err1" "${RCHECK_DIR}/00check.log"
    grep -n "$err2" "${RCHECK_DIR}/00check.log"
    grep -n "$err3" "${RCHECK_DIR}/00check.log"
    grep -n "$err4" "${RCHECK_DIR}/00check.log"
    grep -n "$err5" "${RCHECK_DIR}/00check.log"
    echo "The package contains errors or warnings. See grep results above."
    printf "\n"
    exit 1
fi
