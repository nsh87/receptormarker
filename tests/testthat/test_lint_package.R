context("Lint package to verify code style")


test_that("code conforms to the style guide", {
  lint_errs <- lintr::lint_package(linters=lintr::with_defaults(
    trailing_whitespace_linter=NULL,
    line_length_linter=lintr::line_length_linter(80)))
  
  lint_err_info <- ""
  trim.leading <- function (x)  sub("^\\s+", "", x)
  if (length(lint_errs) > 0) {
    # Some linters failed, build a variable to hold the linter error info:
    k <- lapply(seq_along(lint_errs), function(y, i) {
        x <- y[[i]]  # Set current list element to x
        linted_file <- trim.leading(x["filename"])
        line_num <- x["line_number"]
        col_num <- x["column_number"]
        err_type <- x["type"]
        lint_msg <- x["message"]
        caret_pos <- as.numeric(x["column_number"]) - 1
        # Construct each lint failure's error message
        err_num <- c("[[", i, "]]")
        err_info <- c(linted_file, ":", line_num, ":", col_num, ": ", err_type,
                      ": ", lint_msg)
        guilty_line <- as.character(x["line"])
        caret_str <- paste0(c(rep(" ", caret_pos), "^", "~"), collapse="")
        paste0(c(err_num, "\n",
                 err_info, "\n",
                 guilty_line, "\n",
                 caret_str, "\n\n")
               )
      }
      , y=lint_errs
    )
    lint_err_info <- lapply(k, function(x) {
      trim.leading(paste0(as.character(unlist(x)), collapse=""))
    })
  }
  expect_that(length(lint_errs), equals(0), info=cat(unlist(lint_err_info)))
})
