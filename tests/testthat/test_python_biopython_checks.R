context("System checks for getting Python and Biopython version numbers")


test_that("assume Python exists: verify version number is sane and not NULL", {
  py_vers <- py_version()
  expect_that(py_vers, is_a("numeric"))
  expect_that(length(py_vers), equals(1))
  expect_that(2 < py_vers && py_vers < 4, is_true())
})


test_that("package can laod the Python and Biopython version from options", {
  expect_equal(getOption("receptormarker.py_version"), py_version())
  expect_equal(getOption("receptormarker.biopy_version"), biopy_version())
})


test_that("the versn() function accurately verifies the application name", {
  cmd <- paste0("python -c \"import platform;",
                "print('Python {}'.format(platform.python_version()))\"",
                collapse=" ")
  versn(cmd, "Python")
  expect_null(versn(cmd, "NotPython"))
})


test_that("checking version of fake application doesn't return error message", {
  cmd <- paste0("python -c \"import madeuplibrary;",
                "print('Python {}'.format(madeuplibrary.__version__)\"",
                collapse=" ")
  expect_null(versn(cmd, "MadeUpLibrary"))
})
