context("Unit test multi_clust and validation for it")

data(fluidigm)
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_cluster_plotting.R' test cases.
load(system.file("extdata", "f_clust.rda", package="receptormarker"))

test_that("making sure argument is acceptable range works properly", {
  arg_list <- list("test", NULL, NA, as.factor(2:10), TRUE, 2, data.frame(a=1),
                   list(1), matrix(3:6))
  lapply(arg_list, function(elem) expect_error(validate_sort_range(elem),
                                               "range of integers"))
  arg_list <- list(-1:10, 0:4, 4:-1, c(3, 4, 1, 5))
  lapply(arg_list, function(elem) expect_error(validate_sort_range(elem),
                                               "greater than one"))
  arg_list <- list(c(2, 2, 3, 4), c(2, 3, 3, 4), c(2, 3, 4, 4), c(2, 3, 2, 4))
  lapply(arg_list, function(elem) expect_error(validate_sort_range(elem),
                                               "duplicate values"))
  expect_identical(validate_sort_range(c(2, 4, 3)), c(2, 3, 4))
  expect_identical(validate_sort_range(c(2, 3, 4)), c(2, 3, 4))
  expect_that(validate_sort_range(2:10), not(throws_error()))
})


test_that("making sure argument is multiClust object works properly", {
  arg_list <- list("test", NULL, NA, as.factor(2:10), TRUE, 2, data.frame(a=1),
                   list(1), matrix(3:6))
  lapply(arg_list, function(elem) expect_error(validate_multi_clust(elem),
                                               "object of class 'multiClust'"))
  expect_error(multi_clust(fluidigm[1:25, ], krange = 2:4, runs = 2),
               "not enough rows of data to evaluate")
  expect_that(validate_multi_clust(f_clust), not(throws_error()))
})


test_that("making sure validation functions work properly in multi_clust", {
  expect_error(multi_clust(NULL), "d")
  expect_error(multi_clust(1), "data.frame or matrix")
  expect_error(multi_clust(fluidigm, krange=1), "range of integers")
  expect_error(multi_clust(fluidigm, runs=-1), "runs")
})


test_that("making sure multiClust object works properly", {
  expect_is(f_clust, "multiClust")
  lapply(2:4, function(i) expect_is(f_clust[["sil"]][[i]], "silhouette"))
  expect_equal(length(f_clust[["sil"]]), 4)
  lapply(2:4, function(i) expect_is(f_clust[["sil_avg"]][[i]], "numeric"))
  expect_equal(length(f_clust[["sil_avg"]]), 4)
  lapply(2:4, function(i) {
    expect_equal(f_clust[["sil_avg"]][[i]],
                 summary(f_clust[["sil"]][[i]])[["avg.width"]])
    })
  lapply(2:4, function(i) expect_is(f_clust[["clust_model"]][[i]],
                                    "kmeans"))
  expect_equal(length(f_clust[["clust_model"]]), 4)
  lapply(2:4, function(i) expect_equal(f_clust[["num_clust"]][[i]], i))
  expect_equal(length(f_clust[["num_clust"]]), 4)
  lapply(2:4, function(i) expect_is(f_clust[["wss"]][[i]], "numeric"))
  expect_equal(length(f_clust[["wss"]]), 4)
  expect_is(f_clust[["clust_gap"]], "clusGap")
  expect_equal(length(f_clust[["clust_gap"]]), 4)
  expect_is(f_clust[["k_best"]], "numeric")
  expect_equal(length(f_clust[["k_best"]]), 1)
})
