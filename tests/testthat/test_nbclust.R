context("Unit test NbClust")

# data prep
data(fluidigm)
f50 <- fluidigm[1:50, ]
data(tcr)
tcr_binary_data <- tcr[3:23]
tcr_binary_data[is.na(tcr_binary_data)] <- 0
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_cluster_plotting.R' test cases. It was generated
# using 'data(fluidigm); f_clust <- fluidigm[1:50, ].
#load(system.file("extdata", "f_clust.rda", package="receptormarker"))
# Load TCR data set with only data points and all data points as 0 or 1.
#load(system.file("extdata", "tcr_binary_data.rda", package="receptormarker"))

test_that("making sure invalid arguments produce errors", {
  expect_error(suppressWarnings(NbClust(f50, max.nc = 5, method = 9)),
               "clustering method")
  expect_error(suppressWarnings(NbClust(max.nc = 5, method = "average")),
               "Data matrix")
  expect_error(suppressWarnings(NbClust(f50, max.nc = 5, method = "average",
                                        distance = 9)),
               "distance")
  expect_error(suppressWarnings(NbClust(f50, max.nc = 5, method = "average",
                                        index = 9)),
               "clustering index")
  expect_error(suppressWarnings(NbClust(f50, method = "average",
                                        min.nc = "no")),
               "non-numeric")
  expect_error(suppressWarnings(NbClust(f50, method = "average",
                                        max.nc = "no")),
               "non-numeric")
  expect_error(suppressWarnings(NbClust(f50, max.nc = 5, method = "average",
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
  expect_that(suppressWarnings(NbClust(f50, max.nc = 5, method = "average")),
              not(throws_error()))
})
