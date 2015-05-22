context("System checks for getting Python and Biopython version numbers")

test_that("Python version number is returned and accurate", {
  versn <- py_version()
  expect_that(versn, is_a('numeric'))
  expect_that(length(versn), equals(1))
  expect_that(2 < versn && versn < 4, is_true())
})

test_that("the versn() function accurately verifies the application name", {
  cmd <- paste0('python -c "import platform;',
                'print(\'Python {}\'.format(platform.python_version()))"',
                collapse=' ')
  versn(cmd, 'Python')
  expect_null(versn(cmd, 'NotPython'))
})

test_that("checking version of fake application doesn't return error message", {
  cmd <- paste0('python -c "import madeuplibrary;',
                'print(\'Python {}\'.format(madeuplibrary.__version__)"',
                collapse=' ')
  expect_null(versn(cmd, 'MadeUpLibrary'))
})