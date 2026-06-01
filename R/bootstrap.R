bootstrap_iteration_2 <- function(ip_object) {
  the_call <- ip_object$the_call
  data <- eval(the_call$data, parent.frame())
  predictions <- ip_object$predictions

  bs_sample <- sample(nrow(data), size = nrow(data), replace = T)
  bs_data <- data[bs_sample, ]
  bs_predictions <- lapply(predictions, function(x) x[bs_sample])

  the_call$data <- bs_data
  the_call$object <- bs_predictions
  the_call$bootstrap <- 0
  the_call$null_model <- FALSE
  the_call$strip_ipt_models <- TRUE

  score <- eval.parent(the_call)$score
  return(score)
}

bootstrap_iteration <- function(data, ip_object) {
  # works by creating a new ipscore object based on the original, where all
  # required data has been resampled (& new ipt & ipc weights)

  bs_sample <- sample(nrow(data), size = nrow(data), replace = T)

  # copy ip_object but resample relevant items
  bs_outcome <- list(
    "observed" = ip_object$outcome$observed[bs_sample],
    "type" = ip_object$outcome$type,
    "time_horizon" = ip_object$outcome$time_horizon
  )
  bs_trt <- list(
    "observed" = ip_object$treatment$observed[bs_sample],
    "propensity_formula" = ip_object$treatment$propensity_formula,
    "treatment_of_interest" = ip_object$treatment$treatment_of_interest,
    "type" = ip_object$treatment$type
  )
  bs_pred <- lapply(ip_object$predictions, function(x) x[bs_sample])

  bs_pseudopop <- list(ids = ip_object$pseudopop$ids[bs_sample])

  # compute iptw on resampled data

  bs_iptw <- get_iptw(
    data = data[bs_sample,],
    score_treatment = bs_trt,
    stable_iptw = grepl("stabilized", ip_object$ipt$method),
    only_weights = TRUE
  )

  # if survival, compute ipcw and survival status at time horizon of sample
  if (ip_object$outcome$type == "survival") {
    bs_ipcw <- get_ipcw(
      cens_formula = ip_object$ipc$cens_formula,
      data = data[bs_sample, ],
      cens_model = ip_object$ipc$method,
      time_horizon = ip_object$outcome$time_horizon,
      only_weights = TRUE
    )
    bs_outcome$status_at_horizon <- ip_object$outcome$status_at_horizon[bs_sample]
  } else {
    bs_ipcw <- NULL
  }

  bs_ip_object <- construct_ip_object(
    outcome = bs_outcome,
    treatment = bs_trt,
    predictions = bs_pred,
    pseudopop = bs_pseudopop,
    ipt = bs_iptw,
    ipc = bs_ipcw,
    metrics = ip_object$metrics
  )

  metrics <- compute_metrics(bs_ip_object)$score
  return(metrics)
}

bootstrap_iteration_long <- function(data_outcome, data_long, ip_object) {
  browser()
  bs_sample <- sample(nrow(data_outcome), size = nrow(data_outcome), replace = T)
  bs_data_outcome <- data_outcome[bs_sample, ]
  sampled_ids <- bs_data_outcome$id

  # assign new ids
  bs_data_outcome$id <- seq_along(sampled_ids)

  bs_data_long1 <- do.call(
    rbind,
    lapply(seq_along(sampled_ids), function(i) {

      tmp <- data_long[data_long$id == sampled_ids[i], ]

      # assign corresponding new bootstrap id
      tmp$id <- i

      tmp
    })
  )

  # very slow above. probably some clever matching/merging?

  # get the outcomes from the bootstrapped sample
  bs_outcome <- list(
    "observed" = ip_object$outcome$observed[bs_sample],
    "type" = ip_object$outcome$type,
    "time_horizon" = ip_object$outcome$time_horizon,
    "status_at_horizon" = ip_object$outcome$status_at_horizon[bs_sample]
  )

  # get IPTW on bootstrapped data
  bs_treatment_long <- extract_treatment(bs_data_long,
                                         ip_object$treatment$propensity_formula,
                                         NA)
  bs_ipt <- get_iptw_long(bs_data_long, bs_treatment_long)

  # get treatment info, which is the 'summarized' info on treatment compliance,
  # i.e. 1 row p patient. Resampled compliant patients stay compliant.
  bs_treatment <- ip_object$treatment
  bs_treatment$observed <- bs_treatment$observed[bs_sample]

  # should be same as ip_object$pseudopop but resampled
  bs_pseudopop <- get_pseudopop(bs_outcome, bs_treatment)

  # get censoring weights of bootstrapped sample
  bs_ipc <- get_ipcw_long(
    cens_formula = ip_object$ipc$cens_formula,
    data_outcome = bs_data_outcome,
    data_long = bs_data_long,
    cens_model = ip_object$ipc$method,
    time_horizon = ip_object$outcome$time_horizon,
    strip_ipt_models = TRUE
  )

  # also resample to predictions
  bs_probabilities <- lapply(ip_object$predictions, function(x) x[bs_sample])


  bs_ip_object <- construct_long_ip_object(
    outcome = bs_outcome,
    treatment = bs_treatment,
    predictions = bs_probabilities,
    pseudopop = bs_pseudopop,
    ipt = bs_ipt,
    ipc = bs_ipc,
    metrics = ip_object$metrics
  )

  metrics <- compute_metrics(bs_ip_object)$score

  return(metrics)

}

# arg data_long only used for longitudinal treatment, in which data is data_outcome
bootstrap <- function(data, ip_object, iterations, progress, data_long, type = 1) {
  if (type == 1) {
    b <- lapply_progress(
      as.list(1:iterations),
      function(x) {
        if (identical(class(ip_object), "ip_score")) {
          return(bootstrap_iteration(data, ip_object))
        } else if ("ip_score_long" %in% class(ip_object)) {
          return(bootstrap_iteration_long(
            data_outcome = data, data_long = data_long, ip_object)
          )
        } else {
          stop("unknown class ", class(ip_object), " found for bootstrapping")
        }
      },
      "bootstrapping",
      progress = progress
    )
  } else {
    b <- lapply_progress(
      as.list(1:iterations),
      function(x) {
        if (identical(class(ip_object), "ip_score")) {
          return(bootstrap_iteration_2(ip_object))
        } else if ("ip_score_long" %in% class(ip_object)) {
          return(bootstrap_iteration_long(
            data_outcome = data, data_long = data_long, ip_object)
          )
        } else {
          stop("unknown class ", class(ip_object), " found for bootstrapping")
        }
      },
      "bootstrapping",
      progress = progress
    )
  }

  # transpose results
  # (iteration > metric > model) -> (metric > model > iteration)

  # for calibration plot:
  # (iteration > metric > [pred/obs, model]) ->
  # (metric > model > iteration > list(pred = , obs = ))
  transposed <- lapply(ip_object$metrics, function(m) {
    P <- lapply(names(ip_object$predictions), function(p) {
      if (m != "calplot") { # 1 numeric result, simple to combine & transpose
        sapply(b, function(i) i[[m]][[p]])
      } else { # calibration plot, consisting of 2 vectors of preds & obs
        lapply(b, function(i) {
          list(
            pred = i[[m]][["pred", p]],
            obs = i[[m]][["obs", p]]
          )
        })
      }
    })
    names(P) <- names(ip_object$predictions)
    P
  })
  names(transposed) <- ip_object$metrics

  # # summarize
  conf.int <- lapply(ip_object$metrics, function(m) {
    CI <- lapply(names(ip_object$predictions), function(p) {
      if (m != "calplot") {
        return(ci(transposed[[m]][[p]], cover = 0.95))
      } else {
        return(NA)
      }
    })
    names(CI) <- names(ip_object$predictions)
    CI
  })
  names(conf.int) <- ip_object$metrics

  list(
    results = conf.int,
    raw = transposed
  )
}


lapply_progress <- function(x, FUN, task_description, progress = TRUE) {
  # same as lapply, but print a progress indicator
  n <- length(x)

  if (progress == FALSE) {
    result <- lapply(as.list(1:n), FUN)
  } else {
    FUN2 <- function(x, i, n) {
      result <- FUN(x)
      cat("\r", task_description, ": ", i, "/", n, "     ")
      return(result)
    }

    result <- lapply(as.list(1:n), function(i) FUN2(x[[i]], i, n))
    cat("\r")
  }

  return(result)
}

ci <- function(values, cover = 0.95) {
  lower <- (1-cover) / 2
  upper <- 1 - lower
  stats::quantile(values, probs = c(lower, upper))
}
