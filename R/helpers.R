is_bootstrap_allowed <- function(bootstrap, iptw, ipcw) {
  # do not allow bootstrap if iptw are given as fixed vector
  if (bootstrap != 0 && !missing(iptw) && !is.function(iptw)) {
    stop("can't bootstrap if iptw are given")
  }

  if (bootstrap != 0 && !missing(ipcw) && !is.function(ipcw)) {
    stop("can't bootstrap if ipcw are given")
  }
}

strip_glm <- function(fit) {
  fit[c(
    "residuals",
    "fitted.values",
    "effects",
    "qr",
    "linear.predictors",
    "weights",
    "prior.weights",
    "data"
  )] <- NULL

  fit
}

strip_cox <- function(fit) {
  fit[c(
    "linear.predictors",
    "residuals",
    "y",
    "model"
  )] <- NULL

  fit
}

is.formula <- function(x) {
  inherits(x, "formula")
}

extract_lhs <- function(data, formula) {

  # set rhs to 1, and lhs to . if its missing.
  formula <- stats::update.formula(formula, . ~ 1)

  lhs <- formula[[2]]

  if (lhs == "." || !is.symbol(lhs))
    stop("Left hand side of treatment formula must be one variable")

  eval(lhs, data, emptyenv())
}

rhs_is_one <- function(formula) {
  # set lhs to 1 if there is none, otherwise supplying ~ 1 will return false
  formula <- stats::update.formula(formula, 1 ~ .)
  identical(formula[[3]], 1)
}


check_missing <- function(arg) {
  is_missing <- eval(call("missing", deparse(substitute(arg))),
                     envir = parent.frame())

  if (is_missing) {
    stop("Argument ", as.name(substitute(arg)), " is missing.", call. = FALSE)
  }
}


check_missing_xor <- function(arg1, arg2) {
  is_missing1 <- eval(call("missing", deparse(substitute(arg1))),
                      envir = parent.frame())
  is_missing2 <- eval(call("missing", deparse(substitute(arg2))),
                      envir = parent.frame())

  if (!xor(is_missing1, is_missing2)) {
    stop(
      "One of arguments ", as.name(substitute(arg1)),
      " and ", as.name(substitute(arg2)),
      " must be given",
      call. = FALSE
    )
  }
}

check_input <- function(arg, class) {
  # check if arg is of type class
  if (!inherits(arg, class)) {
    stop(
      sprintf("%s must be of class '%s'.", class),
      call. = FALSE
    )
  }

}



simulate_time_to_event <- function(n, constant_baseline_haz, LP) {
  u <- stats::runif(n)
  -log(u) / (constant_baseline_haz * exp(LP))
}


make_named_list <- function(object, substituteobject) {
  # this function converts a list like list(my_model, model2)
  # to a named list list("my_model" = my_model, "model2" = model2)
  # or an object my_model to list("my_model" = my_model)
  # this gives the models recognizable names

  # exprs <- as.list(substituteobject)
  expr_to_name <- function(x) {
    if (is.character(x)) {
      x
    } else {
      paste(deparse(x, width.cutoff = 20, nlines = 1), collapse = " ")
    }
  }

  # user typed object = list(model1, model2) in function arguments
  if (is.call(substituteobject) &&
      identical(substituteobject[[1]], as.name("list"))) {

    exprs <- as.list(substituteobject)[-1]

    if (is.null(names(object))) {
      newnames <- sapply(exprs, expr_to_name)
    } else {
      newnames <- names(object)

      for (i in seq_along(newnames)) {
        if (newnames[i] == "") {
          newnames[i] <- expr_to_name(exprs[[i]])
        }
      }

    }
    names(object) <- newnames
    return(object)
  }

  # user typed object = x, where x is a list previously defined
  if (identical(class(object), "list")) {
    n <- names(object)
    if (is.null(n)) {
      names(object) <- paste0("model.", 1:length(object))
    } else {
      for (i in seq_along(object)) {
        if (n[i] == "") {
          names(object)[i] <- paste0("model.", i)
        }
      }
    }
    return(object)
  }

  # user typed object = x, where x is a single model
  object <- list(object)
  names(object) <- expr_to_name(substituteobject)
  return(object)
}
