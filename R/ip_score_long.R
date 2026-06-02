# assuming:
# - temperal ordering: P, L0, A0, Y1, L1, A1, Y2, ...
# - visit 0 is on time 0
# - visit times are same for everyone
# - visits are not necessarily equi-spaced?
# - treatment before time 0 is 0.
# - data_long has column 'id' indicating subject id's
# - data_long had column 'visit_id' indicating visit id, starting from 0


#' @export
wide_to_long <- function(df, baseline_variables = c("id"), wide_variables, visit_times,
                         outcome_times) {
  stopifnot("must have column 'id' for subject ids in df" = "id" %in% names(df))
  stopifnot("id must be unique in df" = nrow(df) == length(unique(df$id)))
  long <- reshape(
    data = df,
    varying = wide_variables,
    v.names = names(wide_variables),
    timevar = "visit_time",
    idvar = "id",
    times = visit_times,
    direction = "long",
    new.row.names = NULL
  )
  for (b in baseline_variables) {
    long[[b]] <- rep(df[[b]], length(visit_times))
  }

  # rearrange outcome times in advance, so we can remove rows with visits after
  # survival time later
  outcome_times_rep <- rep(outcome_times, length(visit_times))
  outcome_times_rep <- outcome_times_rep[order(long$id, long$visit_time)]

  cols <- c(baseline_variables, "visit_time", names(wide_variables))
  long <- long[order(long$id, long$visit_time), cols]

  long <- long[long$visit_time <= outcome_times_rep, ]
  rownames(long) <- NULL
  long
}

lag_vec <- function(x, n = 1, fill = 0) {
  nx <- length(x)
  if (n >= nx) {
    return(rep(fill, nx))
  }
  c(rep(fill, n), head(x, -n))
}

#' @export
add_lag_terms <- function(df, var, lag = 1, fill = 0) {
  for (l in lag) {
    lag_name <- paste0(var, "_lag_", l)
    df[[lag_name]] <- unsplit(
      lapply(
        split(df[[var]], df$id),
        lag_vec, n = l, fill = fill
      ),
      df$id
    )
  }
  return(df)
}



#' @param strip_ipt_models If set to TRUE (default), the models for the IPT
#' and IPC-weights are stripped of unnecessary data. Set to FALSE if you plan
#' to do extensive diagnostics on the fitted IPT/IPC models. The resulting
#' `ip_score` object will use quite a lot more memory.
#' @export
ip_score_long <- function(probabilities, data_outcome, data_long,
                          treatment_formula, treatment_of_interest,
                          metrics = c("auc", "brier", "oeratio", "calplot"),
                          visit_times, time_horizon, cens_model = "KM",
                          cens_formula = ~ 1, null_model = TRUE,
                          bootstrap = 0, bootstrap_progress = TRUE,
                          iptw, ipcw, quiet = FALSE,
                          strip_ipt_models = TRUE) {

  # checking inputs...

  check_missing(probabilities)
  check_missing(data_outcome)
  check_missing(data_long)
  check_missing(treatment_formula)
  check_missing(treatment_of_interest)
  check_missing(visit_times)
  check_missing(time_horizon)

  stopifnot(
    "data_outcome requires columns id, time, and status " =
      all(c("id", "time", "status") %in% names(data_outcome))
  )
  stopifnot(
    "data_long should have columns id and visit_time at least" =
      all(c("id", "visit_time") %in% names(data_long))
  )

  if (cens_model == "KM")
    stopifnot("censoring model must be ~ 1 if modeled via KM" =
                rhs_is_one(cens_formula))

  cens_model <- match.arg(cens_model, choices = c("cox", "KM"))

  # do not allow bootstrap if iptw are given as fixed vector
  is_bootstrap_allowed(bootstrap, iptw, ipcw)

  stopifnot(length(unique(data_outcome$id)) == nrow(data_outcome))
  stopifnot(length(unique(data_long$id)) == nrow(data_outcome))
  stopifnot(!is.unsorted(data_long$id)) # triple negation... id must be sorted
  stopifnot(!is.unsorted(data_outcome$id))

  stopifnot(setequal(data_long$visit_time, visit_times))


  # assert:
  # - visits may not be skipped (unless censored). I.e. all visits before
  #   survtime should be in data.

  n_visits <- length(visit_times)
  stopifnot("visit times must be ordered chronologically" =
              all(order(visit_times) == 1:n_visits))
  stopifnot("first visit time must be at t=0" = visit_times[0] == 0)
  stopifnot("last visit time must be before time horizon" =
              max(visit_times) < time_horizon)



  max_visit_time <- tapply(data_long$visit_time, data_long$id, max)
  stopifnot(
    "there are subjects with visit time beyond their survival time" =
      all(max_visit_time <=
            data_outcome$time[match(names(max_visit_time), data_outcome$id)])
  )


  if (is.character(treatment_of_interest) && treatment_of_interest == "always")
    treatment_of_interest <- rep(1, n_visits)
  if (is.character(treatment_of_interest) && treatment_of_interest == "never")
    treatment_of_interest <- rep(0, n_visits)

  score_outcome <- extract_outcome(data_outcome,
                                   substitute(survival::Surv(time, status)),
                                   time_horizon = time_horizon)

  # compute IPT weights. First get a temporary score_treatment_long, which
  # is not used after this.
  score_treatment_long <- extract_treatment(data_long, treatment_formula, NA)
  score_ipt <- get_iptw_long(data_long, score_treatment_long,
                             treatment_of_interest, visit_times,
                             strip_ipt_models, iptw, data_outcome)

  # build a dataframe with 1 row per subject, and whether they are compliant
  # to longitudinal treatment strategy or not. Having early event/censor
  # does not make you non-compliant as long as you are compliant before survtime
  data_flat <- data.frame(
    id = unique(data_long$id),
    ipt = score_ipt$weights,
    trt = faithful_to_trt(data_long, treatment_formula, treatment_of_interest,
                          visit_times)
  )
  names(data_flat)[3] <- as.character(treatment_formula[[2]])

  # we simplified the problem from longitudinal treatment to single treatment
  # (compliant or not). We can reuse a lot of code from ip_score() this way.
  score_treatment <- extract_treatment(data_flat, treatment_formula, TRUE)
  score_pseudopop <- get_pseudopop(score_outcome, score_treatment)

  probabilities <- make_named_list(probabilities, substitute(probabilities))
  score_predictions <- get_predictions(probabilities, data_flat)

  # censoring weights is a bit more tedious, requiring the long data again.
  score_ipc <- get_ipcw_long(cens_formula, data_outcome, data_long,
                       cens_model, time_horizon, visit_times, strip_ipt_models,
                       ipcw)

  if (null_model) {
    score_predictions <- fit_null(score_pseudopop, score_outcome,
                                  score_predictions, score_ipt, score_ipc)
  }

  # Combines all the required elements into 1 object
  ip_object <- construct_long_ip_object(
    outcome = score_outcome,
    treatment = score_treatment,
    predictions = score_predictions,
    pseudopop = score_pseudopop,
    ipt = score_ipt,
    ipc = score_ipc,
    metrics = metrics
  )

  # use the compute_metrics function to compute metrics. All longitudinal
  # information is 'simplified' to a point treatment (compliance/non-compliance)
  # and time dependent weights are combined into 1 weight.
  # Thus, we do not require a special compute_metrics function for long trts.
  ip_object <- compute_metrics(ip_object)

  # do bootstrap
  if (bootstrap > 0) {
    matchcall <- match.call()
    call_env <- parent.frame()

    bs <- bootstrap(ip_object, matchcall, call_env, bootstrap, bootstrap_progress)
    ip_object <- add_to_ip_object(ip_object, "bootstrap", bs, after = 1)
    ip_object <- add_to_ip_object(ip_object, "bootstrap_iterations", bootstrap)
  }

  # trt_of_interest was TRUE, where treatment options was compliant or not
  # compliant in the 'flat' data.frame. Set it back to the long variant as
  # specified for more informative printing
  ip_object$treatment$treatment_of_interest <- treatment_of_interest

  ip_object <- add_to_ip_object(ip_object, "quiet", quiet)

  return(ip_object)
}

faithful_to_trt <- function(data_long, treatment_formula,
                            treatment_of_interest, visit_times) {

  trt_var <- as.character(treatment_formula[[2]])
  visit_id <- match(data_long$visit_time, visit_times)
  as.vector(tapply(
    data_long[[trt_var]] == treatment_of_interest[visit_id],
    data_long$id,
    FUN = function(x) all(x, na.rm = TRUE)
  ))
}

construct_long_ip_object <- function(outcome, treatment, predictions, pseudopop,
                                     ipt, ipc, metrics) {
  ip_object <- construct_ip_object(
    outcome = outcome,
    treatment = treatment,
    predictions = predictions,
    pseudopop = pseudopop,
    ipt = ipt,
    ipc = ipc,
    metrics = metrics
  )
  class(ip_object) <- c("ip_score_long", "ip_score")
  ip_object
}

threshold_weights <- function(weights, quantile_bound) {
  upper_bound <- quantile(weights, quantile_bound)[[1]]
  weights[weights > upper_bound] <- upper_bound
  return(weights)
}

get_iptw_long <- function(data_long, score_treatment, treatment_of_interest,
                          visit_times, strip_model = TRUE, iptw, data_outcome) {

  # if user specified weights themselves:
  if (!missing(iptw)) {
    ipt <- list()
    manualiptw <- handle_specified_ip_long(iptw, data_outcome, data_long)
    ipt$method <- manualiptw$method
    ipt$weights <- manualiptw$ipw
    return(ipt)
  }
  # if not:

  ipt_visit <- get_iptw(data_long, score_treatment, stable_iptw = FALSE,
                        only_weights = FALSE, strip_model = strip_model)

  # set NA when deviating from trt of interest, and set 1 when trt of interest
  # is NA.
  trt_var <- as.character(score_treatment$treatment_column)
  visit_id <- match(data_long$visit_time, visit_times)
  required_trt <- treatment_of_interest[visit_id]
  actual_trt <- data_long[[trt_var]]

  ipt_visit$weights <- ifelse(
    is.na(required_trt),
    1,
    ifelse(
      required_trt == actual_trt,
      ipt_visit$weights,
      NA
    )
  )

  # the above line computes the visit IPT weights. We require the patient
  # ITP weights to be stored in the $weights part.
  ipt_visit$visit_weights <- ipt_visit$weights
  ipt_visit$weights <- as.vector(tapply(
    ipt_visit$visit_weights, data_long$id, FUN = prod
  ))

  return(ipt_visit)
}

get_ipcw_long <- function(cens_formula, data_outcome, data_long,
                          cens_model, time_horizon, visit_times,
                          strip_ipt_models = TRUE, ipcw) {

  # if user specified weights themselves:
  if (!missing(ipcw)) {
    ipc <- list()
    manualipcw <- handle_specified_ip_long(ipcw, data_outcome, data_long)
    ipc$method <- manualipcw$method
    ipc$weights <- manualipcw$ipw
    return(ipc)
  }

  if (cens_model == "KM") {
    # for KM, we can use data_outcome
    ipc <- get_ipcw(Surv(time, status) ~ 1, data_outcome,
                    cens_model = "KM", time_horizon = time_horizon)
  } else if (cens_model == "cox") {
    # for cox, it could be time dependent variables, so we split
    survintervals <- survival::survSplit(
      formula = Surv(time, status) ~ .,
      data = data_outcome,
      cut = visit_times)
    stopifnot(nrow(survintervals) == nrow(data_long))
    data_combined <- cbind(survintervals, data_long)

    if ("ipscore_longcensored" %in% names(data_combined)) {
      stop("Please do not use ipscore_longcensored as one of the columns. This
           name is reserved for internal use.")
    }
    # set censor indicator = 1 in last row of censored patients
    data_combined$ipscore_longcensored <- with(
      data_combined,
      status == 0 & !duplicated(id, fromLast = TRUE))

    full_cens_formula <- stats::update.formula(
      old = cens_formula,
      Surv(tstart, time, ipscore_longcensored) ~ .
    )

    cens_model <- survival::coxph(full_cens_formula, data_combined, model = TRUE)

    bh <- survival::basehaz(cens_model, centered = FALSE)
    cumhaz.fun <- stats::stepfun(bh$time, c(0, bh$hazard))

    # compute the cumulative hazard contribution between each visit
    cumhaz_start <- cumhaz.fun(data_combined$tstart)
    cumhaz_end <- cumhaz.fun(pmin(data_combined$time, time_horizon))
    lp <- stats::predict(cens_model, newdata = data_combined, type = "lp")
    contribution <- -(cumhaz_end - cumhaz_start) * exp(lp)

    # take the sum and compute probability of uncensored
    cumhaz <- tapply(contribution, data_combined$id, FUN = sum)
    prob_uncensor <- exp(cumhaz)

    weight <- ifelse(
      data_outcome[, "status"] == FALSE & data_outcome[, "time"] < time_horizon,
      0,
      1 / prob_uncensor)


    ipc <- list()
    ipc$method <- "cox"
    ipc$cens_formula <- full_cens_formula
    if (strip_ipt_models) {
      ipc$model <- strip_cox(cens_model)
    } else {
      ipc$model <- cens_model
    }
    ipc$weights <- weight

  } else {
    print("censoring model not implemented")
  }
  return(ipc)
}


handle_specified_ip_long <- function(iptcw, data_outcome, data_long, type = "iptw") {
  # if the user specified something in the iptw/ipcw argument, e.g. a vector of
  # weights or a function:
  if (is.vector(iptcw, mode = "numeric")) {

    method <- "weights manually specified"
    ipw <- iptcw
  } else if (is.function(iptcw)) {

    method <- "weights specified via function"

    ipw <- do.call(iptcw, args = list("data_outcome" = data_outcome,
                                      "data_long" = data_long))

    if ( ! (is.vector(ipw, mode = "numeric") && length(ipw) == nrow(data_outcome))) {

      stop("function specified in ", type, " did not return a numeric vector ",
           "of length data_outcome")
    }
  } else {

    stop("argument ", type, " must be missing, a numeric vector of weights, ",
         "or a function.")

  }

  return(list("method" = method, "ipw" = ipw))

}
