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
    if (is.null(a) || !(a %in% c(TRUE, FALSE)) || !(class(a) == "logical")) {
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


#' @title Validate that the arg is either a numeric data.frame or matrix
#' @description An internal function that raises an error if the argument is not
#' either a \emph{data.frame} or \emph{matrix}. Also, If all columns are not
#' numeric, it will raise an error.
#' @param d A \emph{data.frame} or \emph{matrix} whose class the function will
#' confirm.
#' @keywords internal
validate_num_data <- function(d) {
  classes <- c("data.frame", "matrix")
  types <- c("numeric", "integer")
  if (!(class(d) %in% classes)) {
    stop("The argument 'd' is not a data.frame or matrix.", call.=FALSE)
  }
  lapply(d,
         function(x) {
           if (!(class(x) %in% types)) {
             stop("The classes of the columns of 'd' are not all numeric.", 
                  call.=FALSE) 
           } else TRUE
         })
  lapply(d,
         function(x) {
           if (identical(unique(x), c(0, 1))) {
             warning("At least one column of 'd' contains only values 0 and 1.",
                     call.=FALSE)
           }
         })
}

#' @title Validate that an argument contains positive integers
#' @description An internal function that raises an error if the argument does
#' not contain positive integers.
#' @param n A named list of items to be checked.
#' @keywords internal
validate_pos_num <- function(n) {
  lapply(1:length(n), 
         function(i) {
           if (!(class(n[[i]]) %in% c("numeric", "integer")) || 
               length(n[[i]]) > 1) {
             err <- paste0("The argument '", names(n)[[i]],
                           "' must be a positive integer", collapse="")
             stop(err, call.=FALSE)
           } else if (n[[i]] < 1) {
             err <- paste0("The argument '", names(n)[[i]],
                           "' must be a positive integer", collapse="")
             stop(err, call.=FALSE)
           } else if (floor(n[[i]]) != n[[i]]) {
             err <- paste0("The argument '", names(n)[[i]],
                           "' must be a positive integer", collapse="")
             stop(err, call.=FALSE)
           }
         })
}
