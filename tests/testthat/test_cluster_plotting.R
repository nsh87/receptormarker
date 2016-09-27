context("Unit test cluster plotting")

data(fluidigm)
# Load pre-computed multiClust object from fluidigm data set. This same object
# should be used in 'test_multi_clustering.R' test cases. It was generated
# using 'data(fluidigm); f_clust <- multi_clust(fluidigm[1:50, ]).
load(system.file("extdata", "f_clust.rda", package="receptormarker"))

test_that("validation functions work properly in wss_plot", {
  expect_error(wss_plot(NULL), "clust_obj")
  expect_error(wss_plot(f_clust, optimal=1), "optimal")
  expect_error(wss_plot(1), "object of class 'multiClust'")
})


test_that("validation functions work properly in clusGap_plot", {
  expect_error(gap_plot(NULL), "clust_obj")
  expect_error(gap_plot(f_clust, optimal=1), "optimal")
  expect_error(gap_plot(1), "object of class 'multiClust'")
})


test_that("validation functions work properly in pca_plot", {
  expect_error(pca_plot(1, f_clust, num_clust=4), "data.frame or matrix")
  expect_error(pca_plot(fluidigm, 1, num_clust=4),
               "object of class 'multiClust'")
  expect_error(pca_plot(fluidigm, f_clust, num_clust=-1), "num_clust")
  expect_silent(pca_plot(fluidigm, f_clust, num_clust = 2:4))
})


test_that("validation functions work properly in sil_plot", {
  expect_error(sil_plot(1, num_clust=4), "object of class 'multiClust'")
  expect_error(sil_plot(f_clust, num_clust=-1), "num_clust")
})


test_that("validation functions work properly in avg_sil_plot", {
  expect_error(avg_sil_plot(NULL), "clust_obj")
  expect_error(avg_sil_plot(f_clust, optimal=1), "optimal")
  expect_error(avg_sil_plot(1), "object of class 'multiClust'")
})


test_that("validation functions work properly in clust_boxplot", {
  expect_error(clust_boxplot(1, f_clust, num_clust=4), "data.frame or matrix")
  expect_error(clust_boxplot(fluidigm, 1, num_clust=4),
               "object of class 'multiClust'")
  expect_error(clust_boxplot(fluidigm, f_clust, num_clust=-1), "num_clust")
})
