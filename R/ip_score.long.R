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


# input:
# - outcome_data (id, time, status)
# - treatment/confounder data, long. Treatment formula needs to make sense here
#   - add_lag function

#' @export
ip_score_long <- function(probabilities, data_outcome, data_long,
                          treatment_formula, treatment_of_interest,
                          metrics = c("auc", "brier", "oeratio", "calplot"),
                          visit_times, time_horizon, cens_model = "KM",
                          cens_formula = ~ 1, null_model = TRUE, quiet = FALSE) {

  # assert:
  # - probabilities and data outcome same length
  # - data outcome has id, time, status
  # - same ids in data outcome and data long
  # - visit times as given should be consistent with visit_times in data_long
  # - visits may not be skipped (unless censored). I.e. all visits before
  #   survtime should be in data.

  # we should make sure that ids in data_outcome and data_long have same ordering
  # here
  n_visits <- length(visit_times)
  stopifnot("visit times must be ordered chronologically" =
              all(order(visit_times) == 1:n_visits))
  stopifnot("first visit time must be at t=0" = visit_times[0] == 0)
  stopifnot("last visit time must be before time horizon" =
              max(visit_times) < time_horizon)


  if (is.character(treatment_of_interest) && treatment_of_interest == "always")
    treatment_of_interest <- rep(1, n_visits)
  if (is.character(treatment_of_interest) && treatment_of_interest == "never")
    treatment_of_interest <- rep(0, n_visits)

  # create get_iptw long fct instead of this
  ipt_visit <- ipt_weights(data_long, treatment_formula)
  ipt_product <- tapply(ipt_visit$weights, data_long$id,
                        FUN = prod)

  data_flat <- data.frame(
    id = unique(data_long$id),
    ipt = ipt_product,
    trt = faithful_to_trt(data_long, treatment_formula, treatment_of_interest,
                          visit_times)
  )
  names(data_flat)[3] <- as.character(treatment_formula[[2]])

  ipt <- get_iptw(iptw = data_flat$ipt)
  # ipt$model <- ipt_visit$model

  outcome <- extract_outcome(data_outcome, substitute(survival::Surv(time, status)),
                             time_horizon = time_horizon)

  treatment <- extract_treatment(data_flat, treatment_formula, TRUE)


  ipc <- get_ipcw_long(cens_formula, data_outcome, data_long,
                       cens_model, time_horizon)

  predictions <- get_predictions(probabilities, data_flat)

  if (null_model) {
    predictions <- fit_null(treatment, outcome, predictions, ipt, ipc)
  }

  ip_object <- construct_long_ip_object(
    outcome = outcome,
    treatment = treatment,
    predictions = predictions,
    ipt = ipt,
    ipc = ipc,
    metrics = metrics
  )
  ip_object <- compute_metrics(ip_object)

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
    FUN = all
  ))
}

construct_long_ip_object <- function(outcome, treatment, predictions, ipt, ipc,
                                     metrics) {
  ip_object <- construct_ip_object(
    outcome = outcome,
    treatment = treatment,
    predictions = predictions,
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

get_ipcw_long <- function(cens_formula, data_outcome, data_long,
                          cens_model, time_horizon) {
  if (cens_model == "KM") {
    # for KM, we can use data_outcome
    ipc <- get_ipcw(Surv(time, status) ~ 1, data_outcome,
                    cens_model = "KM", time_horizon = time_horizon)
  } else if (cens_model == "cox") {
    survintervals <- survival::survSplit(
      formula = Surv(time, status) ~ .,
      data = data_outcome,
      cut = 0:4)
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
    # ipc$model <- cens_model
    ipc$weights <- weight

  } else {
    print("censoring model not implemented")
  }
  return(ipc)
}
