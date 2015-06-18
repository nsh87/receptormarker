context("Test internal functions that validate other function arguments")

test_that("making sure arguments are TRUE/FALSE works properly", {
  arg_list <- list(fake_param=TRUE, another_param=23)
  expect_error(validate_true_false(arg_list), "another_param")
  expect_error(validate_true_false(list(fake_param=NULL), "fake_param"))
  expect_error(validate_true_false(list(fake_param=NA), "fake_param"))
  arg_list <- list(fake_param=TRUE, another_param=FALSE)
  expect_that(validate_true_false(arg_list), not(throws_error()))
})


test_that("making sure arguments are not NULL works properly", {
  arg_list <- list(fake_param=NULL, another_param=42)
  expect_error(validate_not_null(arg_list), "fake_param")
  arg_list <- list(fake_param=TRUE, another_param=23, third_param=c(1:10),
                   data_frame=c(5:15), fifth_param="test characters",
                   sixth_param=NA)
  expect_that(validate_not_null(arg_list), not(throws_error()))
})
