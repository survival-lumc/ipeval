# pretty print, respecting output width, adding \n at the end
pp <- function(...) {
  txt <- paste0(c(...), collapse = "")
  cat(paste0(strwrap(txt), "\n"))
}


#' @export
print.ip_score <- function(x, ...) {
  if (x$quiet != TRUE) {
    assumptions(x)
  }
  numeric_metrics <- x$metrics[!grepl("^calplot", x$metrics)]

  if (!is.null(x$bootstrap)) {
    for (metric in numeric_metrics) {
      cat("\n", metric, "\n\n", sep = "")
      tab <- data.frame(model = names(x$predictions))
      tab[[metric]] <- x$score[[metric]]
      tab$lower <- sapply(x$bootstrap$results[[metric]], function(x) x[[1]])
      tab$upper <- sapply(x$bootstrap$results[[metric]], function(x) x[[2]])
      print(tab, digits = 3, row.names = FALSE)
    }
  } else {
    cat("\n")
    tab <- data.frame(model = names(x$predictions))
    for (metric in numeric_metrics) {
      tab[[metric]] <- x$score[[metric]]
    }
    print(tab, digits = 3, row.names = FALSE)
  }

  if (any(grepl("calplot", x$metrics))) {
    plot(x, ...)
  }
}

#' @export
plot.ip_score <- function(x, ...) {
  # this plotting function should ideally be more customizable,
  # i.e. show/hide legend, colors, xlim, ylim, ....

  models <- names(x$predictions)

  plot(1, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       asp = 1,
       xlab = "Predicted", ylab = "Observed")
  graphics::title(
    main = paste0("Calibration had everyone followed treatment ",
                  x$treatment$treatment_of_interest),
    col.sub = "#404040",
    cex.sub = 0.8
  )
  graphics::abline(0, 1, col = "black")
  colors <- grDevices::adjustcolor(
    rep(grDevices::palette()[-1], length.out = length(models)),
    alpha.f = 0.8
  )

  for (i in seq_along(models)) {
    graphics::lines(
      x = x$score$calplot[["pred", models[i]]],
      y = x$score$calplot[["obs", models[i]]],
      type = "o",
      col = colors[i],
      lw = 2
    )
  }
  graphics::legend("topleft",
         legend = models,
         col    = colors,
         lty    = 1,
         lwd    = 1,
         pch    = 1,
         bty    = "n")
  if (!is.null(x$bootstrap)) {
    for (m in models) {
      plot(1, type = "n",
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Predicted", ylab = "Observed",
           asp = 1)
      graphics::title(
        main = paste0("Calibration plot for ", m),
        sub = paste0("Calibration had everyone followed treatment ",
                     x$treatment$treatment_of_interest),
        col.sub = "#404040",
        cex.sub = 0.8
      )
      for (i in 1:x$bootstrap_iterations) {
        graphics::lines(
          x = x$bootstrap$raw$calplot[[m]][[i]]$pred,
          y = x$bootstrap$raw$calplot[[m]][[i]]$obs,
          type = "o",
          col = "darkgrey"
        )
      }
      graphics::abline(0, 1, col = "black")
      graphics::lines(
        x = x$score$calplot[["pred", m]],
        y = x$score$calplot[["obs", m]],
        type = "o",
        col = "blue", lw = 2,
      )
      graphics::legend("topleft",
             legend = c("bootstrap iteration", "original (CF) data"),
             col    = c("darkgrey", "blue"),
             lty    = 1,
             lwd    = c(1,2),
             pch    = c(1,1),
             bty    = "n")
    }
  }

}


assumptions <- function(x) {

  if (!is.null(x$ipt$confounders)) {
    confounders <- paste0("{", paste0(x$ipt$confounders, collapse = ", "), "}")
  } else {
    confounders <- "{unknown}"
  }

    if (!is.null(x$ipc$cens_formula)) {
    cv <- all.vars(x$ipc$cens_formula)[-c(1,2)]
    cens_vars <- paste0("{", paste0(cv, collapse = ", "), "}")
  } else {
    cens_vars <- "{unknown}"
  }

  trt_col <- x$treatment$treatment_column
  trt_val <- x$treatment$treatment_of_interest
  pp <- x$pseudopop$ids
  horizon <- x$outcome$time_horizon


  pp("Estimation of the performance of the prediction model in a pseudopopulation
     where everyone's treatment ", trt_col, " was set to ", trt_val, ".")

  if (x$outcome$type == "binary") {
    pp("The pseudopopulation ($pseudopop) is constructed from ", sum(pp),
       " (", round(mean(pp)*100,1),
       "%) subjects who indeed received treatment level ", trt_val,
       ". These subjects are reweighted to represent the full target population
       under a hypothetical intervention in which everyone received this
       treatment level.")
  } else {
    pp("The pseudopopulation ($pseudopop) is constructed from ", sum(pp),
       " (", round(mean(pp)*100,1),
       "%) subjects who indeed received treatment level ", trt_val,
       " and remained uncensored till time=", horizon,
       ". These subjects are reweighted to represent the full target population
       under a hypothetical intervention in which everyone received this
       treatment level and remained uncensored till time=", horizon, ".")
  }

  pp("The following assumptions must be satisfied for correct inference:")

  pp("")
  pp("Causal assumptions:")
  pp("")

  pp("- Conditional exchangeability: after adjustment for the covariates used to construct
  the inverse probability of treatment weights (IPTW), i.e., ", confounders,
  ", there are no unmeasured confounders of treatment assignment and outcome.")

  pp("- Conditional positivity: the probability of receiving treatment level ",
     trt_val, " should be greated than zero for each combination of the
     variables ", confounders, " that is observed in the full population. The
     distribution of IPT-weights can be assessed with
     $ipt$weights[$pseudopop$ids].")

  pp("- Consistency: the observed outcome under the received treatment equals
     the potential outcome under that treatment. This includes the assumption of
     no interference between subjects.")

  if (x$outcome$type == "survival") {
    if (x$ipc$method == "KM") {
      pp("- Independent censoring. The censoring mechanism is completely
         independent of the outcome process.")
      pp("- Positivity for censoring: requires that the probability of remaining
         uncensored till time=", horizon, " is greater than zero. The
         distribution of IPC-weights can be assessed with
         $ipc$weights[$pseudopop$ids].")
    } else if (x$ipc$method == "cox") {
      pp("- Conditionally independent censoring: conditional on variables ",
         cens_vars, ", censoring is independent of the outcome process.")
      pp("- Conditional positivity for censoring: requires that for all observed
         combinations of the covariate variables ", cens_vars, " the probability
         of remaining uncensored till time=", horizon, " is greater than zero.
         The distribution of IPC-weights can be assessed with
         $ipc$weights[$pseudopop$ids].")
    } else {
      pp("- The censoring mechanism is adjusted for by the manually specified
         inverse probability of censoring weights (IPCW).")
    }
  }

  pp("")
  pp("Modeling assumptions:")
  pp("")

  if (grepl("binomial", x$ipt$method)) {

    pp("- Correctly specified propensity model. Estimated treatment model is")
    pp(print_model(x$ipt$model), ". See also $ipt$model.")

  } else if (grepl("multinomial", x$ipt$method)) {

    pp("- Correctly specified propensity model. The multinomial treatment model
       is estimated using the nnet package. See also $ipt$model.")

  } else {

    pp("- The supplied inverse probability of treatment weights (IPTW) are
       assumed to be valid.")

  }

  if (grepl("stabilized", x$ipt$method)) {

    pp("* Stabilized weights were used. See also $ipt$stable_model.
      Pseudopopulation weights ($ipt$weights) are the probability of
      treatment from $ipt$stable_model divided by probability of treatment
         from $ipt$model.")
  }

  if (x$outcome$type == "survival") {

    switch(x$ipc$method,

    KM = {
      pp("- The censoring distribution was estimated nonparametrically using
             the Kaplan-Meier estimator. The probability of remaining uncensored
            is ")
      pp("  ", print_km(x$ipc$model, horizon),
         ". See also $ipc$model.")
    },
    cox = {
      pp("- Correctly specified censoring model. The estimated censoring
      model is ")
      pp(" ", print_censor_model(x$ipc$model), ". See also $ipc$model.")
    },

    pp("- The supplied inverse probability of censoring weights (IPCW) are
       assumed to be valid.")
    )
  }

  pp("")
  pp("Performance estimates:")
}


print_model <- function(model) {
  link <- model$family$link
  lhs_var <- model$formula[[2]]

  lhs <- paste0(link, "(", lhs_var, ")")

  var_names <- names(model$coefficients)

  coef <- round(unname(model$coefficients), 2)

  rhs <- paste(coef, var_names, sep = "*", collapse = " + ")
  rhs <- gsub("*(Intercept)", "", rhs, fixed = TRUE)
  rhs <- gsub("+ -", "- ", rhs, fixed = TRUE)

  formula <- paste(lhs, rhs, sep = " = ")
  formula
}

print_censor_model <- function(cox) {
  if (is.null(cox$coefficients)) {
    LP <- "0"
  } else {
    var_names <- names(cox$coefficients)
    coef <- round(unname(cox$coefficients), 2)
    LP <- paste(coef, var_names, sep = "*", collapse = " + ")
  }
  paste0("P(C > t) = C_0(t)^exp(", LP, ")")
}

print_km <- function(km, time_horizon) {
  paste0("P(C > ", time_horizon, ") = ",
         round(stats::predict(km, times = time_horizon), 2))
}
