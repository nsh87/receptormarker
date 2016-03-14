context("Unit test multi_clust and validation for it")

data(iris)
data(fluidigm)
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_cluster_plotting.R' test cases. It was generated
# using 'data(fluidigm); f_clust <- multi_clust(fluidigm[1:50, ]).
load(system.file("extdata", "f_clust.rda", package="receptormarker"))
# Load TCR data set with only data points and all data points as 0 or 1.
load(system.file("extdata", "tcr_binary_data.rda", package="receptormarker"))
# Create a contrived boolean data set
contrived_bool <- data.frame(matrix(0, nrow=90, ncol=9))
contrived_bool[1:20, 1] <- 1
contrived_bool[21:35, 9] <- 1
contrived_bool[36:90, 5] <- 1

test_that("argument validation of the multiClust object works properly", {
  arg_list <- list("test", NULL, NA, as.factor(2:10), TRUE, 2, data.frame(a=1),
                   list(1), matrix(3:6))
  lapply(arg_list, function(elem) expect_error(validate_multi_clust(elem),
                                               "object of class 'multiClust'"))
  expect_that(validate_multi_clust(f_clust), not(throws_error()))
})


test_that("argument validation functions work properly in multi_clust", {
  expect_error(multi_clust(NULL), "d")
  expect_error(multi_clust(1), "data.frame or matrix")
  expect_error(multi_clust(fluidigm, krange=1), "range of integers")
  expect_error(multi_clust(fluidigm, runs=-1), "runs")
})


test_that("you indeed have a multiClust object for testing with", {
  expect_is(f_clust, "multiClust")
})

test_that("the 'sil' slot in the pre-generated multiClust validates", {
  expect_equal(length(f_clust@sil), 15)
})

test_that("the 'sil_avg' slot in the pre-generated multiClust validates", {
  expect_equal(length(f_clust@sil_avg), 15)
})

test_that("the 'clust_model' slot in the pre-generated multiClust validates", {
  expect_is(f_clust@clust_model[[2]], "kmeans")
  expect_equal(length(f_clust@clust_model), 15)
})
  
test_that("the 'num_clust' slot in the pre-generated multiClust validates", {
  lapply(2:4, function(i) expect_equal(f_clust@num_clust[[i]], i))
  expect_equal(length(f_clust@num_clust), 15)
})

test_that("the 'wss' slot in the pre-generated multiClust validates", {
  expect_equal(length(f_clust@wss), 15)
})

test_that("the 'clust_gap' slot in the pre-generated multiClust validates", {
  expect_equal(length(f_clust@clust_gap), 4)
})

test_that("the 'k_best' slot in the pre-generated multiClust validates", {
  expect_equal(length(f_clust@k_best), 1)
})

test_that("multiClust object can be generated with boolean data", {
  expect_that(multi_clust(tcr_binary_data), not(throws_error()))
})

test_that("multiClust object can be generated with non-boolean data", {
  expect_that(multi_clust(fluidigm[1:40, ]), not(throws_error()))
})

# The data sets used below are tested directly with NbClust in test_nbclust.R,
# however multi_clust() is supposed to detect boolean data itself and pass the
# correct arguments to NbClust; this basically is an end-to-end test to ensure
# that happens correctly.
test_that("multi_clust() correctly handles contrived boolean data", {
  k_best <- multi_clust(contrived_bool, krange = 2:7)@k_best
  expect_identical(k_best, 3)
})

test_that("multi_clust() correctly handles tcr boolean data", {
  k_best <- multi_clust(tcr_binary_data, krange = 2:7)@k_best
  expect_identical(k_best, 2)
})

test_that("multi_clust() correctly handles non-boolean iris data", {
  k_best <- multi_clust(iris[, 1:4], krange = 2:7)@k_best
  expect_identical(k_best, 3)
})

test_that("multi_clust() correctly handles non-boolean fluidigm data", {
  k_best <- multi_clust(fluidigm[1:50, ], krange = 2:7)@k_best
  expect_identical(k_best, 3)
})
