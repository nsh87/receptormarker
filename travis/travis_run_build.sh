#!/usr/bin/env sh

# Build and check the package --as-cran
echo "Building with R CMD build"
R CMD build .
export PKG_TARBALL=$(Rscript -e 'pkg <- devtools::as.package("."); cat(paste0(pkg$package, "_", pkg$version, ".tar.gz"));')

# Check the package
echo "Testing with: R CMD check ${PKG_TARBALL} --as-cran"
R CMD check "${PKG_TARBALL}" --as-cran
# Check for warnings
export RCHECK_DIR=$(Rscript -e 'cat(paste0(devtools::as.package(".")$package, ".Rcheck"))')
grep -q -R "WARNING" "${RCHECK_DIR}/00check.log"; RETVAL=$?
