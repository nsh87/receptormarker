[![Travis-CI Build Status](https://travis-ci.org/nsh87/receptormarker.svg?branch=master)](https://travis-ci.org/nsh87/receptormarker)

## Installation

While this package is in development it can be installed directly from GitHub
using the `devtools` package:

```R
devtools::install_github('nsh87/receptormarker')
```

Once version 1.0 is complete the package will be submitted to CRAN.

## Development

Fork the repository and submit changes with a pull request.

To get started, make sure you have the latest version of R and RStudio
installed. It is strongly suggested that you have installed the following
packages:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
```

Load **receptormarker.Rproj** to open the project in RStudio, then load the
code:

```R
devtools::load_all()
```

Be sure to follow this [code style](http://r-pkgs.had.co.nz/r.html#style "Hadley Wickham's Modified Google R Style Guide")
 and [document your functions](http://r-pkgs.had.co.nz/man.html "Documenting
with Roxygen2").  If you create any functions that are internal and not of
interest to most users, be sure to document them with `@keywords internal` to
exclude them from the package index. You should explicitly
[define functions to export](http://r-pkgs.had.co.nz/namespace.html#exports "Namespacing in R")
to NAMESPACE with Roxygen2's `@external` tag. Build .Rd files from Roxygen2 comments with:

```R
devtools::document()
```

It is strongly suggested that you [automate
testing](http://r-pkgs.had.co.nz/tests.html "Writing Tests for R") of all
functions written. Run all tests with:

```R
devtools::test()
```

Before submitting a pull request you should build the documentation, run tests
and check that the package builds. This can be done with a single command:

```R
devtools::check()
```
