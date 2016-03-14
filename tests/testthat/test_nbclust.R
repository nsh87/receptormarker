context("Unit test NbClust")

# data prep
data(fluidigm)
f40 <- fluidigm[1:40, ]
data(iris)
rep_row <- function(d, n) {
  new_d <- c()
  for (i in 1:n) {
    new_d <- rbind(new_d, d)
  }
  return(new_d)
}
clust_data <- rbind(rep_row(c(1, rep(0, 8)), 20),
                    rep_row(c(rep(0, 8), 1), 15),
                    rep_row(c(rep(0, 4), 1, rep(0, 4)), 55))
rownames(clust_data) <- 1:90
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_cluster_plotting.R' test cases. It was generated
# using 'data(fluidigm); f_clust <- fluidigm[1:50, ].
load(system.file("extdata", "f_clust.rda", package="receptormarker"))
# Load TCR data set with only data points and all data points as 0 or 1.
load(system.file("extdata", "tcr_binary_data.rda", package="receptormarker"))

test_that("making sure invalid arguments produce errors", {
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

test_that("making sure NbClust object can be generated with boolean data", {
  expect_that(suppressWarnings(NbClust(tcr_binary_data,
                               distance = "binary",
                               max.nc = 5,
                               method = "average")),
              not(throws_error()))
})

test_that("making sure NbClust object can be generated with non-boolean data", {
  expect_that(suppressWarnings(NbClust(f40, max.nc = 5, method = "average")),
              not(throws_error()))
})

test_that("making sure NbClust object picks right number of clusters", {
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
  nb_best <- suppressWarnings(NbClust(iris[1:4],
                                      min.nc = 3,
                                      index = "alllong",
                                      max.nc = 7,
                                      method = "average"))
  best <- aggregate(nb_best[["Best.nc"]][1, ], 
                    by = list(nb_best[["Best.nc"]][1, ]), length)
  index <- which.max(best[[2]])
  k_best <- best[index, 1]
  expect_identical(k_best, 3)
  nb_best <- suppressWarnings(NbClust(clust_data,
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
