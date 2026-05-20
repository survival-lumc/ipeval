flip_surv_event <- function(formula) {
  lhs <- formula[[2]]

  lhs[[3]] <- call("!", lhs[[3]])   # negate event indicator
  formula[[2]] <- lhs
  formula
}

ipc_weights <- function(data, formula, type, time_horizon) {
  # perhaps use riskRegression::ipcw()
  if (type == "KM")
    stopifnot(rhs_is_one(formula))

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)

  p_uncensored <- switch(
    type,
    KM = {
      fit <- prodlim::prodlim(formula, data = data, reverse = TRUE)
      p_not_censor <- stats::stepfun(fit$time, c(1, fit$surv), right = TRUE)
      list(
        model = fit,
        probability = p_not_censor(pmin(y[, "time"], time_horizon))
      )
    },
    cox = {
      # coxph has no reverse argument, need to flip it manually
      # not sure whether this handles event/censor ties correctly
      flipped_form <- flip_surv_event(formula)

      fit <- survival::coxph(flipped_form, data = data, model = TRUE, x = TRUE)
      list(
        model = fit,
        probability = 1-predict_cox(fit, data, pmin(y[, "time"], time_horizon))
      )
    },
    stop("cens.model ", type, " not implemented")
  )

  # if censored before time horizon, weight is NA,
  # else, weight is 1/probability uncensored at event/time horizon
  w <- ifelse(
    y[, "status"] == 0 & y[, "time"] < time_horizon,
    NA,
    1 / p_uncensored$probability
  )

  list(
    model = p_uncensored$model,
    weights = unname(w)
  )
}
