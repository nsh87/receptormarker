context("Unit test NbClust")

# Data prep
data(iris)
data(fluidigm)
f40 <- fluidigm[1:40, ]
contrived_bool <- data.frame(matrix(0, nrow=90, ncol=9))
contrived_bool[1:20, 1] <- 1
contrived_bool[21:35, 9] <- 1
contrived_bool[36:90, 5] <- 1
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_cluster_plotting.R' test cases. It was generated
# using 'data(fluidigm); f_clust <- multi_clust(fluidigm[1:50, ]).
load(system.file("extdata", "f_clust.rda", package="receptormarker"))
# Load TCR data set with only data points and all data points as 0 or 1.
load(system.file("extdata", "tcr_binary_data.rda", package="receptormarker"))

test_that("invalid arguments to NbClust produce errors", {
  expect_error(suppressWarnings(NbClust(f40, max.nc = 5, method = 9)),
               "clustering method")
  expect_error(suppressWarnings(NbClust(max.nc = 5, method = "average")),
               "Data matrix")
  expect_error(suppressWarnings(NbClust(f40, max.nc = 5, method = "average",
                                        distance = 9)),
               "distance")
  expect_error(suppressWarnings(NbClust(f40, max.nc = 5, method = "average",
                                        index = 9)),
               "clustering index")
  expect_error(suppressWarnings(NbClust(f40, method = "average",
                                        min.nc = "no")),
               "non-numeric")
  expect_error(suppressWarnings(NbClust(f40, method = "average",
                                        max.nc = "no")),
               "non-numeric")
  expect_error(suppressWarnings(NbClust(f40, max.nc = 5, method = "average",
                                        min.nc = 7)),
               "difference between the minimum")
})

test_that("NbClust object can be generated with boolean data", {
  expect_that(suppressWarnings(NbClust(tcr_binary_data,
                               distance = "binary",
                               max.nc = 5,
                               method = "average")),
              not(throws_error()))
})

test_that("NbClust object can be generated with non-boolean data", {
  expect_that(suppressWarnings(NbClust(f40, max.nc = 5, method = "average")),
              not(throws_error()))
})

test_that("NbClust object picks right number of clusters with fluidigm", {
  nb_best <- suppressWarnings(NbClust(f40,
                                      min.nc = 2,
                                      index = "alllong",
                                      max.nc = 7,
                                      method = "average"))
  best <- aggregate(nb_best[["Best.nc"]][1, ], 
                    by = list(nb_best[["Best.nc"]][1, ]), length)
  index <- which.max(best[[2]])
  k_best <- best[index, 1]
  expect_identical(k_best, 3)
})

test_that("NbClust object picks right number of clusters with tcr binary", {
  nb_best <- suppressWarnings(NbClust(tcr_binary_data,
                                      min.nc = 2,
                                      index = "alllong",
                                      max.nc = 7,
                                      distance = "binary",
                                      method = "average"))
  best <- aggregate(nb_best[["Best.nc"]][1, ], 
                    by = list(nb_best[["Best.nc"]][1, ]), length)
  index <- which.max(best[[2]])
  k_best <- best[index, 1]
  expect_identical(k_best, 2)
})

test_that("NbClust object picks right number of clusters with iris", {
  nb_best <- suppressWarnings(NbClust(iris[, 1:4],
                                      min.nc = 3,
                                      index = "alllong",
                                      max.nc = 7,
                                      method = "average"))
  best <- aggregate(nb_best[["Best.nc"]][1, ], 
                    by = list(nb_best[["Best.nc"]][1, ]), length)
  index <- which.max(best[[2]])
  k_best <- best[index, 1]
  expect_identical(k_best, 3)
})
  
test_that("NbClust object picks right number of clusters with contrived bool", {
  nb_best <- suppressWarnings(NbClust(contrived_bool,
                                      min.nc = 2,
                                      index = "alllong",
                                      max.nc = 7,
                                      distance = "binary",
                                      method = "average"))
  best <- aggregate(nb_best[["Best.nc"]][1, ], 
                    by = list(nb_best[["Best.nc"]][1, ]), length)
  index <- which.max(best[[2]])
  k_best <- best[index, 1]
  expect_identical(k_best, 3)
})
