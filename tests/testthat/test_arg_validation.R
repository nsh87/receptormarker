context("Test internal functions that validate other function arguments")

test_that("validation of TRUE/FALSE works properly", {
  arg_list <- list(fake_param=TRUE, another_param=23)
  expect_error(validate_true_false(arg_list), "another_param")
  arg_list <- list(fake_param=TRUE, another_param=1)
  expect_error(validate_true_false(arg_list), "another_param")
  arg_list <- list(fake_param=TRUE, another_param=0)
  expect_error(validate_true_false(arg_list), "another_param")
  expect_error(validate_true_false(list(fake_param=NULL), "fake_param"))
  expect_error(validate_true_false(list(fake_param=NA), "fake_param"))
  arg_list <- list(fake_param=TRUE, another_param=FALSE)
  expect_silent(validate_true_false(arg_list))
})


test_that("not NULL argument validation works properly", {
  arg_list <- list(fake_param=NULL, another_param=42)
  expect_error(validate_not_null(arg_list), "fake_param")
  arg_list <- list(fake_param=TRUE, another_param=23, third_param=c(1:10),
                   data_frame=data.frame(x=1, y=1:10, let="abc"),
                   fifth_param="test characters", sixth_param=NA)
  expect_silent(validate_not_null(arg_list))
})


test_that("argument validation for multi_clust works properly", {
  arg_list <- list("test", NULL, NA, as.factor(2:10), TRUE, 2, data.frame(a=1),
                   list(1), matrix(3:6))
  lapply(arg_list, function(elem) expect_error(validate_k_range(elem),
                                               "range of integers"))
  arg_list <- list(-1:10, 0:4, 4:-1, c(3, 4, 1, 5))
  lapply(arg_list, function(elem) expect_error(validate_k_range(elem),
                                               "greater than one"))
  arg_list <- list(c(2, 2, 3, 4), c(2, 3, 3, 4), c(2, 3, 4, 4), c(2, 3, 2, 4))
  lapply(arg_list, function(elem) expect_error(validate_k_range(elem),
                                               "duplicate values"))
  expect_silent(validate_k_range(2:10))
})
