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
