# assuming:
# - temperal ordering: P, L0, A0, Y1, L1, A1, Y2, ...
# - visit 0 is on time 0
# - visit times are same for everyone
# - visits are not necessarily equi-spaced?
# - treatment before time 0 is 0.
# - data_long has column 'id' indicating subject id's
# - data_long had column 'visit_id' indicating visit id, starting from 0


#' @export
wide_to_long <- function(df, baseline_variables, wide_variables, visit_times,
                         outcome_times) {
  # we sometimes want a time dependent variable to be available in all rows too
  both_i <- baseline_variables %in% unlist(wide_variables)
  both_vars <- baseline_variables[both_i]
  baseline_variables <- baseline_variables[!both_i]

  long <- reshape(
    data = df,
    varying = wide_variables,
    v.names = names(wide_variables),
    timevar = "visit_time",
    idvar = baseline_variables,
    times = visit_times,
    direction = "long",
    new.row.names = NULL
  )
  for (b in both_vars) {
    long[[b]] <- rep(df[[b]], length(visit_times))
  }

  # rearrange outcome times in advance, so we can remove rows with visits after
  # survival time later
  outcome_times_rep <- rep(outcome_times, length(visit_times))
  outcome_times_rep <- outcome_times_rep[order(long$id, long$visit_time)]

  cols <- c(baseline_variables, both_vars, "visit_time", names(wide_variables))
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
ip_score_long <- function(probabilities, data_outcome, data_long, visit_times,
                          time_horizon, treatment_formula, treatment_of_interest,
                          metrics = c("auc", "brier", "oeratio", "calplot"),
                          null_model = TRUE) {

  # assert:
  # - probabilities and data outcome same length
  # - data outcome has id, time, status
  # - same ids in data outcome and data long
  # - visit times should be correct


  # we should make sure that ids in data_outcome and data_long have same ordering
  # here
  stopifnot(visit_times[0] == 0)
  stopifnot(max(visit_times) < time_horizon)
  n_visits <- length(visit_times)

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
    trt = faithful_to_trt(data_long, treatment_of_interest)
  )
  names(data_flat)[3] <- as.character(treatment_formula[[2]])

  outcome <- extract_outcome(data_outcome, substitute(survival::Surv(time, status)),
                             time_horizon = time_horizon)

  treatment <- extract_treatment(data_flat, treatment_formula, TRUE)

  ipt <- get_iptw(iptw = data_flat$ipt)
  ipt$model <- ipt_visit$model
  ipc <- get_ipcw(ipcw = rep(1, nrow(data_flat)))

  # ipt$weights <- threshold_weights(ipt$weights, 0.99)

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

  ip_object <- add_to_ip_object(ip_object, "quiet", FALSE)
  ip_object <- add_to_ip_object(ip_object, "treatment_formula", ipt_visit$model)

  # TODO: print in assumptions # patients that satisfy treatment strategy
  # also in point trt!

  return(ip_object)
}

faithful_to_trt <- function(data_long, treatment_of_interest) {
  # TODO: replace A by trt var, lhs from treatment formula
  # TODO: "visit_time" not hardcoded?
  # TODO: data_long[["visit"]] + 1 doesnt work when visit times are not 0,1,2,..
  tapply(
    data_long$A == treatment_of_interest[data_long[["visit_time"]] + 1],
    data_long$id,
    FUN = all
  )
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
