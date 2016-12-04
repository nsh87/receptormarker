context("Test the convergence pipeline")

# Data prep
data(tcr)
tcr_small <- tcr[0:10, ]


test_that("results of convergence pipeline are parsed correctly", {
  c <- convergence(tcr_small, "seqs")
  expect_s4_class(c, "convergenceGroups")
  # Check that the data frame of results comes back correctly
  expect_equal(nrow(c@groups), 6)
  expect_equal(names(c@groups), c("num_items", "group_name", "X1"))
  # Check that the data frame of network info comes back correctly
  expect_equal(nrow(c@network), 6)
  expect_equal(names(c@network), c("node1", "node2", "type"))
  # Check that there are no rows with empty values in groups
  c@groups[c@groups == ""] <- NA
  expect_identical(anyNA(c@groups), FALSE)
  # Check that there are no rows with empty values in network
  c@network[c@network == ""] <- NA
  expect_identical(anyNA(c@network), FALSE)
  # Make sure we have sequences in various columns
  expect_true(is.element("CASSEAGGQDYGNEQFF", c@groups[, "X1"]))
  expect_true(is.element("CRG-CASSEAGGQDYGNEQFF", c@groups[, "group_name"]))
  expect_true(is.element("CASPGGWTGGGNEQFF", c@network[, "node1"]))
  expect_true(is.element("CASNTLGAGGREQYF", c@network[, "node2"]))
  expect_true(is.element("singleton", c@network[, "type"]))
})
