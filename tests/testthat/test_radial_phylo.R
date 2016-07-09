context("Unit test radial phylogram")

data(tcr)
tcr <- tcr[1:25, ]
tcr_evil <- tcr; tcr_evil[3, "seqs"] <- "EVIL123"
tcr_onecol <- data.frame(tcr[, "seqs"]); names(tcr_onecol) <- "seqs"
tcr_vect <- as.character(tcr_onecol[, "seqs"])


test_that("validation of user data and seqence column works in all cases", {
  # Data validation
  expect_error(radial_phylo(d=TRUE),
               "'d' must be a data.frame or a character vector")
  expect_error(radial_phylo(d=NA),
               "'d' must be a data.frame or a character vector")
  expect_error(radial_phylo(d=c(1:10)),
               "'d' must be a data.frame or a character vector")
  expect_error(radial_phylo(d=NULL),
               "cannot be NULL")
  expect_error(radial_phylo(d=c("CASSYPQT", "CA1235679")),
               "Sequences must only contain characters from A-Z")
  expect_error(radial_phylo(d=tcr_evil, "seqs"),
               "Sequences must only contain characters from A-Z")

  # Sequence column
  expect_error(radial_phylo(tcr_vect, seqs_col="seqs"),
               "'seqs_col' expected to be NULL")
  expect_error(radial_phylo(tcr_vect, seqs_col=1),
               "'seqs_col' expected to be NULL")
  expect_error(radial_phylo(tcr_vect, seqs_col=NA),
               "'seqs_col' expected to be NULL")
  expect_error(radial_phylo(tcr, seqs_col="col_not_exists"),
               "Data does not have column named 'col_not_exists'")
  expect_error(radial_phylo(tcr, seqs_col=NA),
               "'seqs_col' must be a column index or a column name")
  expect_error(radial_phylo(tcr, seqs_col=NULL),
               "'seqs_col' must be a column index or a column name")
  expect_error(radial_phylo(tcr, seqs_col=c(1, 2)),
               "'seqs_col' must represent a single column")
  expect_error(radial_phylo(tcr, seqs_col=82),
               "'seqs_col' is greater than the number of columns of data")
  expect_error(radial_phylo(tcr_onecol, seqs_col="col_not_exists"),
               "Data does not have column named 'col_not_exists'")
  expect_error(radial_phylo(tcr_onecol, seqs_col=NA),
               "'seqs_col' must be a column index or a column name")
  expect_error(radial_phylo(tcr_onecol, seqs_col=NULL),
               "'seqs_col' must be a column index or a column name")
  expect_error(radial_phylo(tcr_onecol, seqs_col=c(1, 2)),
               "'seqs_col' must represent a single column")
  expect_error(radial_phylo(tcr_onecol, seqs_col=82),
               "'seqs_col' is greater than the number of columns of data")
})

  
test_that("validation of 'font_size' works", {
  expect_error(radial_phylo(tcr, "seqs", font_size=1000),
               "'font_size' must be an integer between 1 and 99")
  expect_error(radial_phylo(tcr, "seqs", font_size=-10),
               "'font_size' must be an integer between 1 and 99")
  expect_error(radial_phylo(tcr, "seqs", font_size="12"),
               "'font_size' must be an integer")
  expect_error(suppressWarnings(radial_phylo(tcr, "seqs", font_size=c(1:10))),
               "'font_size' must be an integer")
  expect_error(radial_phylo(tcr, "seqs", font_size=NA),
               "'font_size' must be a real value")
  expect_error(radial_phylo(tcr, "seqs", font_size=TRUE),
               "'font_size' must be an integer")
})


test_that("validation of 'canvas_size' works", {
  expect_error(radial_phylo(tcr, "seqs", canvas_size="1000"),
               "'canvas_size' is invalid")
  expect_error(radial_phylo(tcr, "seqs", canvas_size=-10),
               "'canvas_size' is invalid")
  expect_error(radial_phylo(tcr, "seqs", canvas_size=NA),
               "'canvas_size' is invalid")
  expect_error(radial_phylo(tcr, "seqs", canvas_size=c(1:10)),
               "'canvas_size' is invalid")
  expect_error(radial_phylo(tcr, "seqs", canvas_size=TRUE),
               "'canvas_size' is invalid")
})

  
test_that("validation of 'rings' works", {
  # Using a data frame
  expect_error(radial_phylo(tcr, "seqs", rings=c("col_fake"=1)),
               "the column name 'col_fake' does not exist in input data")
  expect_error(radial_phylo(tcr, "seqs", rings=c(1:3)),
               "'rings' must be a named vector, or NULL")
  expect_error(radial_phylo(tcr, "seqs", rings=TRUE),
               "'rings' must be a named vector, or NULL")
  expect_error(radial_phylo(tcr, "seqs", rings=NA),
               "'rings' must be a named vector, or NULL")
  expect_error(radial_phylo(tcr, "seqs", rings="test"),
               "elements in the argument 'rings' must be named")
  expect_error(radial_phylo(tcr, "seqs", rings=42),
               "elements in the argument 'rings' must be named")

  # Using a single-column data frame (not allowed with 'rings')
  expect_error(radial_phylo(tcr_onecol, "seqs", rings=c("col_doesnt_exist"=1)),
               "'d' must have more than one column")
  expect_error(radial_phylo(tcr_onecol, "seqs", rings=c(seqs="CASYTS")),
               "'d' must have more than one column")

  # Using a character vector (not allowed)
  expect_error(radial_phylo(tcr_vect, rings=c("col_doesnt_exist"=1)),
               "'d' expected to be a data.frame")
  expect_error(radial_phylo(tcr_vect, rings=c(seqs="CASYTS")),
               "'d' expected to be a data.frame")
})


radial_phylo_widget <- c("radial_phylo", "htmlwidget")  # Class of the widget
test_that("radial phylograms work using R", {
  # Vanilla phylograms
  r <- radial_phylo(tcr, "seqs", fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr, 2, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_onecol, 1, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_onecol, "seqs", fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_vect, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))

  # With condense option
  r <- radial_phylo(tcr, "seqs", condense=TRUE, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr, 2, condense=TRUE, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_onecol, 1, condense=TRUE, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_onecol, "seqs", condense=TRUE, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr_vect, condense=TRUE, fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))

  # With condense and rings
  r <- radial_phylo(tcr, "seqs", condense=TRUE, rings=c(FOXP3=1), fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr, 2, condense=TRUE, rings=c(FOXP3=1), fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))

  # With just rings
  r <- radial_phylo(tcr, "seqs", rings=c(FOXP3=1), fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
  r <- radial_phylo(tcr, 2, rings=c(FOXP3=1), fast=TRUE)
  expect_that(r, is_a(radial_phylo_widget))
})


test_that("radial phylograms work using Biopython", {
  biopy_existence <- check_bio_python(level="warn")
  if(!(is.null(biopy_existence)) && biopy_existence == TRUE) {
    # Vanilla phylograms
    r <- radial_phylo(tcr, "seqs")
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr, 2)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_onecol, 1)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_onecol, "seqs")
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_vect)
    expect_that(r, is_a(radial_phylo_widget))

    # With condense option
    r <- radial_phylo(tcr, "seqs", condense=TRUE)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr, 2, condense=TRUE)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_onecol, 1, condense=TRUE)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_onecol, "seqs", condense=TRUE)
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr_vect, condense=TRUE)
    expect_that(r, is_a(radial_phylo_widget))

    # With condense and rings
    r <- radial_phylo(tcr, "seqs", condense=TRUE, rings=c(FOXP3=1))
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr, 2, condense=TRUE, rings=c(FOXP3=1))
    expect_that(r, is_a(radial_phylo_widget))

    # With just rings
    r <- radial_phylo(tcr, "seqs", rings=c(FOXP3=1))
    expect_that(r, is_a(radial_phylo_widget))
    r <- radial_phylo(tcr, 2, rings=c(FOXP3=1))
    expect_that(r, is_a(radial_phylo_widget))
  } else {
    expect_that(radial_phylo(tcr, "seqs"),
                gives_warning("install Python's Biopython"))
  }
})
