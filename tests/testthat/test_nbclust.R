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

test_that("making sure argument is acceptable range works properly", {
  
})

test_that("multiClust object can be generated with boolean data", {
  expect_that(NbClust(tcr_binary_data,
                      distance = "binary",
                      max.nc = 5,
                      method = "average"),
              not(throws_error()))
})

test_that("multiClust object can be generated with non-boolean data", {
  expect_that(multi_clust(fluidigm[1:40, ]), not(throws_error()))
})
