#' @title Validate than an object is of class \emph{convergenceGroups}
#' @description An internal function that raises an error if the argument is not
#' of class \code{\link{convergenceGroups-class}} or if the \code{groups} are
#' empty (i.e. if there are no convergence groups). 
#' @param convergence_obj An item to be checked for class membership.
#' @keywords internal
validate_convergence_clust <- function(convergence_obj) {
  if (class(convergence_obj) != "convergenceGroups") {
    err <- paste0(c("The argument 'convergence_obj' must be an object of",
                    "class 'convergenceGroups'"),
                  sep=" ")
    stop(err, call.=FALSE)
  }
  # Get the data.frame of convergence groups and check if there are 0 rows
  if (nrow(convergence_obj@groups) < 1) {
    err <- paste0(c("There are no convergence groups in",
                    "'convergence_obj@groups'"),
                  sep=" ")
    stop(err, call.=FALSE)
  }
}


#' @title Validate that an object is of class \code{\link{multiClust-class}}
#' @description An internal function that raises an error if the argument is not
#' of class \code{\link{multiClust-class}}.
#' @param clust_obj An item to be checked for class membership.
#' @keywords internal
validate_multi_clust <- function(clust_obj) {
  if (class(clust_obj) != "multiClust") {
    stop("The argument 'clust_obj' must be an object of class 'multiClust'",
         call.=FALSE)
  }
}


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


#' @title Validate protein or DNA sequences
#' @description An internal function that creates an error if sequences contain
#' characters outisde the alphabet.
#' @param seqs A character vector of sequences.
#' @keywords internal
validate_sequences <- function(seqs) {
  # Make sure sequences are only alpha characters
  seqs_col_err <- "Sequences must only contain characters from A-Z and a-z"
  g <- grepl("[^A-Za-z]", as.character(seqs))
  if (sum(g) > 0) {
    stop(seqs_col_err, call.=FALSE)
  }
}

#' @title Validate that the arg is either a numeric data.frame or matrix
#' @description An internal function that raises an error if the argument is not
#' either a \emph{data.frame} or \emph{matrix}. Also, if all columns are not
#' numeric, it will raise an error.
#' @param d A \emph{data.frame} or \emph{matrix} whose class (and the classes
#' of each of its columns) the function will confirm.
#' @keywords internal
validate_num_data <- function(d) {
  classes <- c("data.frame", "matrix")  # Acceptable classes for 'd'
  types <- c("numeric", "integer")  # Acceptable types for each column of 'd'
  if (!(class(d) %in% classes)) {
    stop("The argument 'd' is not a data.frame or matrix.", call.=FALSE)
  }
  lapply(d, function(x) {
    if (!(class(x) %in% types)) {
      stop("The classes of the columns of 'd' are not all numeric.", 
           call.=FALSE) 
    }
  })
  for (i in ncol(d)) {
    col <- d[, i]
    uniq <- unique(col)
    if (all(uniq %in% 0:1)) {
      message("At least one column of 'd' contains only values 0 and 1.")
      break
    }
  }
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

#' @title Validate than an argument is a single, positive integer
#' @description An internal function that raises an error if the argument is not
#' a positive integer.
#' @param n An item to be checked.
#' @keywords internal
validate_single_pos_num <- function(n) {
  if (is.na(n)) {
    stop("Argument 'n' must be a real value", call.=FALSE)
  } else if (class(n) != "numeric") {
    stop("Argument 'n' must be an integer", call.=FALSE)
  } else if (length(n) != 1) {
    stop("Argument 'n' must be a single integer", call.=FALSE)
  } else if(n == 0) {
    stop("Argument 'n' must be greater than 0", call.=FALSE)
  } else if (grep("^[0-9]+$", n) != 1) {
    stop("Argument 'n' must be an integer greater than 0",
         call.=FALSE)
  }
}

#' @title Validate the krange parameter of the multi_clust function
#' @description An internal function that raises an error if the argument is not
#' a positive, non-duplicate range of integers beginning > 1. It is sorted into
#' ascending order.
#' @param n_range An item to be checked to make sure it is a valid range.
#' @keywords internal
validate_k_range <- function(range) {
  if (!(class(range) %in% c("integer", "numeric")) || length(range) <= 1) {
    stop("The argument 'krange' must be a range of integers.", call.=FALSE)
  } else if (any(range <= 1)) {
    stop("The argument 'krange' must contain only values greater than one.", 
         call.=FALSE)
  } else if (anyDuplicated(range) != 0) {
    stop("The argument 'krange' must not contain any duplicate values.",
         call.=FALSE)
  }
}
