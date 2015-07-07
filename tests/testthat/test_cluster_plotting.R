context("Unit test cluster plotting")

tryCatch({
  set.seed(1)
  data(fluidigm)
  fluidigm <- fluidigm[1:25,]
  f_clust <- multi_clust(fluidigm, krange = 2:4)
},
finally = {
  set.seed(NULL)  # Turn off seed
}
)

test_that("making sure validation functions work properly in wss_plot", {
  expect_error(wss_plot(NULL), "clust_obj")
  expect_error(wss_plot(f_clust, optimal=1), "optimal")
  expect_error(wss_plot(1), "object of class 'multiClust'")
})


test_that("making sure wss plots properly", {
  hash <- "22341ee84dd9450dd59dabec185c137f" # used optimal = FALSE
  hash1 <- "6f19d768eb28da7688fcfa24a21bba0c" # used optimal = TRUE
  file_path <- tempfile(pattern="wss_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  wss_plot(f_clust, optimal=FALSE)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  file_path <- tempfile(pattern="wss_plot", fileext=".png")
  png(file_path)
  wss_plot(f_clust, optimal=TRUE)
  dev.off()
  test_hash1 <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
  expect_equal(hash1, test_hash1)
})


test_that("making sure validation functions work properly in clusGap_plot", {
  expect_error(clusGap_plot(NULL), "clust_obj")
  expect_error(clusGap_plot(f_clust, optimal=1), "optimal")
  expect_error(clusGap_plot(1), "object of class 'multiClust'")
})


test_that("making sure clusGap plots properly", {
  hash <- "661a488755796b7e8f8d2da7d68241e7" # used optimal = FALSE
  hash1 <- "9ab1de7304fb78e8c816a22bab5f8d4a" # used optimal = TRUE
  file_path <- tempfile(pattern="clusGap_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  clusGap_plot(f_clust, optimal=FALSE)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  file_path <- tempfile(pattern="clusGap_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  clusGap_plot(f_clust, optimal=TRUE)
  dev.off()
  test_hash1 <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
  expect_equal(hash1, test_hash1)
})


test_that("making sure validation functions work properly in pca_plot", {
  expect_error(pca_plot(1, f_clust, num_clust=4), "data.frame or matrix")
  expect_error(pca_plot(fluidigm, 1, num_clust=4),
               "object of class 'multiClust'")
  expect_error(pca_plot(fluidigm, f_clust, num_clust=-1), "num_clust")
})


test_that("making sure pca plots properly", {
  hash <- "a4b9a2bc0549b9086bfa8e93f9159eca"
  file_path <- tempfile(pattern="pca_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  pca_plot(fluidigm, f_clust, num_clust=4)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
})


test_that("making sure validation functions work properly in sil_plot", {
  expect_error(sil_plot(1, num_clust=4), "object of class 'multiClust'")
  expect_error(sil_plot(f_clust, num_clust=-1), "num_clust")
})


test_that("making sure silhouette plots properly", {
  hash <- "307a3e18adb81469e9be84c320cc5cf1"
  file_path <- tempfile(pattern="sil_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  sil_plot(f_clust, num_clust=4)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
})


test_that("making sure validation functions work properly in avg_sil_plot", {
  expect_error(avg_sil_plot(NULL), "clust_obj")
  expect_error(avg_sil_plot(f_clust, optimal=1), "optimal")
  expect_error(avg_sil_plot(1), "object of class 'multiClust'")
})


test_that("making sure avg_sil plots properly", {
  hash <- "2dbc56d557beb6b0b92dc9b111fbdd50" # used optimal = FALSE
  hash1 <- "04b2dc3bef38f8bc90e5a41386849092" # used optimal = TRUE
  file_path <- tempfile(pattern="avg_sil_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  avg_sil_plot(f_clust, optimal=FALSE)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  file_path <- tempfile(pattern="avg_sil_plot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  avg_sil_plot(f_clust, optimal=TRUE)
  dev.off()
  test_hash1 <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
  expect_equal(hash1, test_hash1)
})


test_that("making sure validation functions work properly in clust_boxplot", {
  expect_error(clust_boxplot(1, f_clust, num_clust=4), "data.frame or matrix")
  expect_error(clust_boxplot(fluidigm, 1, num_clust=4),
               "object of class 'multiClust'")
  expect_error(clust_boxplot(fluidigm, f_clust, num_clust=-1), "num_clust")
})


test_that("making sure clust_boxplots properly", {
  hash <- "f40043e4a25bdbcc5c650531f23c74be"
  file_path <- tempfile(pattern="clust_boxplot", fileext=".png")
  png(file_path, width=480, height=480, antialias="none")
  clust_boxplot(fluidigm, f_clust, num_clust=4)
  dev.off()
  test_hash <- digest::digest(file.info(file_path)[["size"]], algo="md5")
  expect_equal(hash, test_hash)
})
