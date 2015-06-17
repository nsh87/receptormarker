#' @title Validate that a list of arguments are \code{TRUE} or \code{FALSE}
#' @description An internal function that raises an error if any of the items in
#' the list is not either \code{TRUE} or \code{FALSE}.
#' @param arg_list A named list containing the variables to check. The names in
#' the list should correspond to the variable names. For example,
#' \code{list(param=param, datum=datum)}.
#' @keywords internal
validate_true_false <- function(arg_list) {
  for (i in c(1:length(arg_list))) {
    a <- arg_list[[i]]
    if (is.null(a) || !(a %in% c(TRUE, FALSE))) {
      err <- paste0("The argument '", names(arg_list)[[i]],
                "' must be TRUE or FALSE", collapse="")
      stop(err, call.=FALSE)
    }
  }
}


#' @title Validate that a list of arguments are not \code{NULL}
#' @description An internal function that raises an error if any of the items in
#' the list is \code{NULL}.
#' @inheritParams validate_true_false
#' @keywords internal
validate_not_null <- function(arg_list) {
  for (i in c(1:length(arg_list))) {
    a <- arg_list[[i]]
    if (is.null(a)) {
      err <- paste0("The argument '", names(arg_list)[[i]],
                    "' cannot be NULL", collapse="")
      stop(err, call.=FALSE)
    }
  }
}
