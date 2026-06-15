#' Convert wide longitudinal data to long format
#'
#' Reshapes repeated measurements stored in wide format into a long-format
#' dataset with one row per subject per visit. Rows corresponding to visits
#' occurring after the subject's outcome time are removed.
#'
#' @param df A data frame containing one row per subject.
#' @param baseline_variables Character vector of baseline variables to retain
#'   and repeat across visits. Must include `"id"`.
#' @param wide_variables Named list mapping long-format variable names to the
#'   corresponding wide-format column names for each visit.
#' @param visit_times Numeric vector giving the visit times corresponding to the
#'   repeated measurements in `wide_variables`.
#' @param outcome_times Numeric vector of outcome or follow-up times, used to
#'   remove visits occurring after a subject's observed follow-up.
#'
#' @returns A data frame in long format containing one row per subject per
#'   observed visit, with columns `id`, `visit_time`, the specified baseline
#'   variables, and the reshaped longitudinal variables.
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#'   id = 1:3,
#'   A0 = c(0, 1, 0),
#'   A1 = c(1, NA, 0),
#'   L0 = c(0.2, -1.1, 0.5),
#'   L1 = c(0.8, NA, 0.1),
#'   time = c(3, 1, 4)
#' )
#'
#' wide_to_long(
#'   df = data,
#'   baseline_variables = "id",
#'   wide_variables = list(
#'     A = c("A0", "A1"),
#'     L = c("L0", "L1")
#'   ),
#'   visit_times = c(0, 2),
#'   outcome_times = data$time
#' )
wide_to_long <- function(df, baseline_variables = c("id"), wide_variables, visit_times,
                         outcome_times) {
  stopifnot("must have column 'id' for subject ids in df" = "id" %in% names(df))
  stopifnot("id must be unique in df" = nrow(df) == length(unique(df$id)))
  long <- stats::reshape(
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
  c(rep(fill, n), utils::head(x, -n))
}

#' Add lagged versions of a longitudinal variable
#'
#' Creates one or more lagged versions of a variable in long-format data. Lagged
#' values are computed separately within each subject and added as new columns.
#' The data must be sorted by id and visit time before calling this function.
#'
#' @param df A long-format data frame containing an `id` column identifying
#'   subjects.
#' @param var Character string giving the name of the variable to lag.
#' @param lag Integer vector specifying the lag(s) to create.
#' @param fill Value used when a lagged observation is unavailable (e.g. at the
#'   first visit). Defaults to 0.
#'
#' @returns The input data frame with additional columns named
#'   `(var)_lag_(lag)` containing the requested lagged values.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   id = c(1, 1, 1, 2, 2, 2),
#'   visit_time = c(0, 1, 2, 0, 1, 2),
#'   A = c(0, 1, 1, 1, 0, 0)
#' )
#'
#' # Add a 1-visit lag of A
#' add_lag_terms(df, "A")
#'
#' # Add both 1- and 2-visit lags
#' add_lag_terms(df, "A", lag = c(1, 2))
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

#'Interventional prediction score for longitudinal treatment strategies
#'
#'Estimates the predictive performance of predictions of binary (or
#'time-to-event) outcomes under longitudinal intervention strategies, by
#'reweighting the data to form a pseudo-population in which every subject was
#'assigned the longitudinal treatment strategy of interest.
#'
#'To form a pseudo-population that represents a setting in which everybody
#'received treatment level 1 during five visits, set `treatment_of_interest` to
#'`c(1,1,1,1,1)`. Any pattern can be chosen, e.g. `c(1,0,1,0,1)` is also valid,
#'as long as the number of visits equals ‘n_visits’. Alternatively, one could
#'also set this argument to "always" or "never" as a shortcut for `rep(1,
#'n_visits)` or `rep(0, n_visits)`. Treatment strategies where treatment is only
#'set at certain visits are also possible via `NA`, e.g. use `c(1,1,NA, NA, NA)`
#'when you want to form a pseudo-population where everybody's treatment levels
#'are set to 1 at the first 2 visits, and their remaining three can be whatever
#'they would have normally been after those first two under the natural course.
#'If treatment is categorical, this should be a character vector denoting the
#'treatment levels of interest, e.g `c(“active”,”control”).
#'
#'The KM censoring distribution is estimated using `prodlim::prodlim(...,
#'reverse = TRUE)`. This correctly estimates the censoring distribution when
#'there are ties between event and censoring times.  When using a Cox model to
#'estimate the censoring distribution, the event indicator is reversed manually. This does
#'not preserve the usual tie-handling convention: in standard survival analysis,
#'censoring is assumed to occur after events at the same time point, but after
#'reversing the indicator the opposite ordering is assumed. A possible
#'workaround is to add a small positive offset (`epsilon`) to all censoring
#'times before fitting the censoring model.
#'
#'Bootstrapping is not possible when manually specifiying the IPTW/IPCW as
#'numeric vectors. If specifying a model user-defined function that computes the
#'ITPW/IPCW given data, it is possible. The given function will be called on
#'each bootstrapped dataset and resulting metrics are used to compute the
#'(empirical) 95\% Cis. More advanced techniques such as thresholding extreme IP
#'weights, can be implemented through a user-defined weight function. The
#'censoring weight returned by this function should be the 1 / probability of
#'remaining uncensored at the time horizon, or at their event time, whichever
#'happens first.
#'
#'@param predictions A numeric vector corresponding to the risk estimates
#'  under evaluation, or a (named) list of multiple numeric vectors for
#'  evaluating and comparing multiple vectors.
#'@param data_outcome A dataframe containing the observed outcomes. It must
#'  consist of 1 row per subject, with columns 'id', 'time', and 'status'. Time
#'  can be continuous, status represents the binary outcome status at ‘time’.
#'@param data_long A dataframe in 'long' format, containing the time-dependent
#'  treatment variable, the (potentially time-dependent) adjustment variables
#'  (confounders) and other time dependent covariates, possibly required for
#'  modeling the censoring mechanism. Baseline covariates can also be included,
#'  in which case they would be repeated for every visit. The data should be in
#'  'long' format, i.e. each subject has one row for each of their visits. Must
#'  include the columns 'id' and 'visit_time'. There should not be any visits
#'  after a subject's follow up time. All subjects should have data at each
#'  visit.
#'@param treatment_formula A formula which indicates the treatment/intervention
#'  (left hand side) and the adjustment variables (right hand side) in
#'  data_long. E.g. A ~ L + A_lag_1. The left hand side can be either a binary
#'  treatment (coded as 0/1 numeric, logical or factor) or a treatment with more
#'  than two categories (coded as factor). The adjustment variables (right hand
#'  side) are used to estimate the inverse probability of treatment weights
#'  (IPTW). All variables on the right hand side must be present in data_long.
#'  The IPTW can also be specified directly using the iptw argument, in which
#'  case the right hand side of this formula is ignored (the left hand side must
#'  still indicate the treatment, i.e. A ~ 1).
#'@param treatment_of_interest A treatment strategy under which the predictions
#'  should be evaluated. Should be in the form of a vector, with one element per
#'  visit. See details.
#'@param metrics A character vector specifying which performance metrics to be
#'  computed. Options are c("auc", "brier", “scaled_brier”, "oeratio",
#'  "calplot").
#'@param visit_times A numeric vector, indicating the times of the visits at
#'  which sequential treatment decisions were made. The first visit must always
#'  be at time 0. Visit times are not necessarily equidistant.
#'@param time_horizon the prediction horizon of interest. Should be after the
#'  last visit.
#'@param cens_model Model for estimating inverse probability of censored weights
#'  (IPCW). Methods currently implemented are Kaplan-Meier ("KM") or Cox
#'  ("cox"), with censoring times derived from the (time,status) variables in
#'  data_outcome, reversing the status indicator, see details. KM is only
#'  supported when the right hand side of cens_formula is 1.
#'@param cens_formula Model formula used for estimating the censoring
#'  probabilities, e.g. ~ x1 + x2. Could consist of only baseline variables but
#'  also of time-dependent variables. All variables specified must be present in
#'  data_long.
#'@param null_model If TRUE a model without covariates (intercept only) is
#'  fitted to data that estimates the same probability for all subjects in data.
#'  The model is fitted using the reweighted data in which all subjects are
#'  'counterfactually' assigned the treatment strategy of interest (using the
#'  IPTW, as estimated using the treatment_formula or as given by the iptw
#'  argument). For censored outcomes, the subjects are also
#'  'counterfactually' uncensored (using the IPCW, as estimated using the
#'  cens_formula, or as given by the ipcw argument). The null model can be used
#'  as a reference (baseline) model.
#'@param bootstrap If this is an integer greater than 0, this indicates the
#'  number of bootstrap iterations, to compute 95\% confidence intervals around
#'  the performance metrics.
#'@param bootstrap_progress if set to TRUE, print a progress bar indicating the
#'  progress of bootstrap procedure.
#'@param iptw A numeric vector with length nrow(data_outcome), containing one
#'  over the probability of being compliant to the treatment strategy of interest
#'  for each subject If iptw is not specified, these
#'  weights are computed using the treatment_formula, but they can be specified
#'  directly via this argument. A user-defined function can also be specified,
#'  which takes as input arguments 'data_outcome' and 'data_long' and returns a
#'  numeric vector of the IPTW weights. See details.
#'@param ipcw A numeric vector with length nrow(data_outcome), containing the
#'  inverse probability of censoring weights at the time horizon, or at their event time, whichever happens
#'  first. If ipcw is not specified, these weights are computed using the
#'  cens_formula, but they can be specified
#'  directly via this argument. A user-defined function can also be specified,
#'  which takes as input arguments 'data_outcome' and 'data_long' and returns a
#'  numeric vector of the IPCW weights. See details.
#'@param quiet  If set to TRUE, don't print assumptions.
#'@param strip_ipt_models If set to TRUE (default), unnecessary components from
#'  the IPT- and IPC-model objects are not stored to save memory. Set to FALSE
#'  if you want to store the full IPT/IPC model objects.
#'
#'@returns An object of class `ip_score`, for which the `print()` and `plot()`
#' methods are implemented. The object is a nested list containing: \itemize{
#'   \item `$score`, contains the estimated predictive performance metrics.
#'   \item `$bootstrap`, if requested, the 95\% confidence intervals of the
#'   performance metrics, and the performance metrics for each individual
#'   bootstrap iteration.
#'   \item `$outcome`, the observed outcomes from data_outcome.
#'   \item `$treatment`, the observed treatment levels from data_long.
#'   \item `$predictions`, the predictions to be evaluated, i.e. the estimated probability
#'   of event under the intervention strategy of interest for each subject.
#'   \item `$ipt`, method, model and inverse probability of treatment weights
#'   (IPTW). These are NA for subjects who are not directly used in the
#'   pseudo-population.
#'   \item `$ipc`, method, model and inverse probability of censoring weights
#'   (IPCW). These are NA for subjects who were censored.
#'   \item `$pseudopop`, binary vector indicating which subjects of the original
#'   population were used to create the pseudo-population, by receiving the
#'   treatment strategy of interest and remaining uncensored, if applicable.
#'   }
#'  The print method summarizes the results and if (quiet = FALSE), prints the
#'  assumptions required for valid inference.
#'
#'@export
#'@seealso \code{\link{ip_score}}
#'
#' Open the corresponding vignette for more extensive examples with
#' \code{vignette("longitudinal", package = "ipeval")}
#'@references Keogh RH, Van Geloven N. Prediction Under Interventions:
#'  Evaluation of Counterfactual Performance Using Longitudinal Observational
#'  Data. Epidemiology. 2024;35(3):329-339.
#'
#'@examples
#' set.seed(5)
#' n <- 1000
#' data <- data.frame(id = 1:n)
#'
#' # 2 visits at t = 0 and t = 2. Time dependent confounding between A and L
#' # Here we generate random survival times between 0 and 6 as an example.
#' # See simulating-data vignette for another example where survival outcomes
#' # depend on time dependent  L and A.
#' data <- within(data, {
#'   L0 <- rnorm(n)
#'   A0 <- rbinom(n, 1, plogis(0.7*L0))
#'   L1 <- rnorm(n, 0.8*L0 - 0.2*A0)
#'   A1 <- rbinom(n, 1, plogis(0.7*L1 + 0.6*A0))
#'   time <- runif(n, 0, 6)
#'   status <- rbinom(n, 1, 0.5)
#' })
#' # A0 is at time t = 0, A1 is at time t = 2
#'
#' # If you have data in this 'wide' form, then usually you can use the convenience
#' # functions to make it suitable for ip_score_long
#' head(data)
#' # note that due to our simulation there may be measurements after a subject was
#' # censored or had event. This will be removed later.
#'
#  # Data preperation for ip_score_long
#' data_outcome <- data[, c("id", "time", "status")]
#' data_long <- wide_to_long(
#'   df = data,
#'   baseline_variables = "id",
#'   wide_variables = list("A" = c("A0", "A1"),
#'                         "L" = c("L0", "L1")),
#'   visit_times = c(0,2),
#'   outcome_times = data$time
#' )
#' # adding treatment level at previous visit time to the data. For the first
#' # visit this is assumed to be 0 (fill)
#' data_long <- add_lag_terms(data_long, var = "A", lag = 1, fill = 0)
#'
#' head(data_long)
#' # note that the measurements that were generated after follow up time have been
#' # removed by wide_to_long().
#'
#' # To get some predictions to evaluate, here is a model that
#' # randomly predicts a number between 0.2 and 0.8
#' random_model <- runif(n, min = 0.2, max = 0.8)
#'
#' # Our treatment at visit i depends on confounder value L at same visit and,
#' # if visit = 2, on treatment level at previous visit. We will illustrate a correctly
#' # specified IPTW model with A ~ L + A_lag_1.
#'
#' # Estimate performance of random predictions in a pseudopopulation in
#' # which everybody was assigned treatment strategy {0, 0}
#' ip_score_long(
#'   predictions = random_model,
#'   data_outcome = data_outcome,
#'   data_long = data_long,
#'   treatment_formula = A ~ L + A_lag_1,
#'   treatment_of_interest = c(0,0),
#'   visit_times = c(0,2),
#'   time_horizon = 4
#' )
#'
#' # Performance of random predictions in a pseudopopulation in
#' # which everybody was assigned treatment level 1 at visit 1, and whatever
#' # they would normally have after that (natural course) at visit 2, and modelling
#' # the censoring distribution with Cox model
#' ip_score_long(
#'   predictions = random_model,
#'   data_outcome = data_outcome,
#'   data_long = data_long,
#'   treatment_formula = A ~ L + A_lag_1,
#'   treatment_of_interest = c(1,NA),
#'   visit_times = c(0,2),
#'   time_horizon = 4,
#'   cens_model = "cox",
#'   cens_formula = ~ A + A_lag_1 + L
#' )
ip_score_long <- function(predictions, data_outcome, data_long,
                          treatment_formula, treatment_of_interest,
                          metrics = c("auc", "brier", "scaled_brier", "oeratio", "calplot"),
                          visit_times, time_horizon, cens_model = "KM",
                          cens_formula = ~ 1, null_model = TRUE,
                          bootstrap = 0, bootstrap_progress = TRUE,
                          iptw, ipcw, quiet = FALSE,
                          strip_ipt_models = TRUE) {
  # checking inputs...

  check_missing(predictions)
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
  stopifnot(!is.unsorted(data_long$id)) # triple negation: id must be sorted
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

  # build a dataframe with 1 row per subject, and whether they are compliant
  # to longitudinal treatment strategy or not. Having early event/censor
  # does not make you non-compliant as long as you are compliant before survtime
  data_flat <- data.frame(
    id = unique(data_long$id),
    trt = faithful_to_trt(data_long, treatment_formula, treatment_of_interest,
                          visit_times)
  )
  names(data_flat)[2] <- as.character(treatment_formula[[2]])

  # we simplified the problem from longitudinal treatment to single treatment
  # (compliant or not). We can reuse a lot of code from ip_score() this way.
  score_treatment <- extract_treatment(data_flat, treatment_formula, TRUE)
  score_pseudopop <- get_pseudopop(score_outcome, score_treatment)

  predictions <- make_named_list(predictions, substitute(predictions))
  score_predictions <- get_predictions(predictions, data_flat)

  # compute IPT/IPC weights. This is done on the longitudinal data. Weights
  # are combined as products to get 1 weight per patient, which represents
  # the weight that a patient is compliant to treatment(in the data_flat sense).
  # First get a temporary score_treatment_long, which is not used after this.
  score_treatment_long <- extract_treatment(data_long, treatment_formula, NA)
  score_ipt <- get_iptw_long(data_long, score_treatment_long,
                             treatment_of_interest, visit_times,
                             strip_ipt_models, iptw, data_outcome)

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
