#' Interventional prediction score
#'
#' Estimates the performance of predictions of binary or time-to-event outcomes under baseline
#' interventions, by reweighting the data to form a pseudo-population in which
#' every subject was assigned the treatment level of interest.
#'
#' When supplying a glm or coxph model as object, the function will try to
#' estimate risks from the model under the treatment level of interest for all
#' subjects in data. If the model does not have the treatment as covariate, it
#' is assumed it already estimates the risk under the treatment level of
#' interest (e.g. because the model was fitted in a population all receiving the treatment level of interest). Alternatively, if the model includes the treatment as covariate,
#' the function estimates the risk under the treatment level of interest for all
#' subjects in data, even if they were assigned an alternative treatment level.
#'
#' All performance metrics are computed on the weighted population mimicking the
#' hypothetical situation where every subject’s treatment level was set to the
#' treatment level of interest (and where nobody was censored before time_horizon). "auc" is area
#' under the (ROC) curve. "brier" is Brier score, ranging from 0 to 1. Scaled
#' brier score is also available (metrics = "scaled_brier"), which expresses the Brier score relative to the null model. For the O/E ratio,
#' the numerator (observed) is the (weighted) fraction of 'observed' events in
#' the pseudo-population, and the denominator (expected) is the (unweighted) mean
#' of risk estimates for all subjects in data. The calplot option generates a
#' calibration plot, with default 8 subgroups. More/less subgroups can be
#' specified by appending “calplot” with a number indicating the number of
#' subgroups, e.g. metrics = "calplot10" for 10 subgroups.
#'
#' The null model estimates are obtained as the weighted mean outcome in the subset observed with the treatment level of interest (and not censored before time_horizon). For time-to-event data, this null prediction could also be
#' computed using a weighted Kaplan-Meier estimator, which would be more
#' efficient, but computationally slower.
#'
#' The censoring distribution is estimated with a Kaplan-Meier estimator implemented using `prodlim::prodlim(...,
#' reverse = TRUE)`. This correctly estimates the censoring distribution when
#' there are ties between event and censoring times.  When using a Cox model to
#' estimate the censoring distribution, the event indicator is reversed. This
#' does not preserve the usual tie-handling convention: in standard survival
#' analysis, censoring is assumed to occur after events at the same time point,
#' but after reversing the indicator the opposite ordering is assumed. A
#' possible workaround is to add a small positive offset (`epsilon`) to all
#' censoring times before fitting the censoring model.
#'
#' Bootstrapping is not possible when manually specifying the IPTW/IPCW as
#' numeric vectors. If specifying a user-defined function that computes the
#' IPTW/IPCW given data, it is possible. The given function will be called on
#' each bootstrapped dataset and resulting metrics are used to compute the 95\%
#' CIs with the percentile method. More advanced techniques, such as thresholding extreme IP weights, can
#' be implemented through user-defined weight function. The censoring weight
#' returned by this function should be the 1 / probability of remaining
#' uncensored till the time_horizon, or till a subject's event time, whichever happens
#' first. For subjects who are censored before the time_horizon, the ipcw can be left at NA.
#'
#' @param object One of the following three options can be used to input the
#'   predictions to be evaluated:
#' \itemize{
#'   \item a numeric vector, corresponding to the risk estimates under evaluation
#'   \item a glm or coxph model, from which the predictions under evaluation can
#'   be derived. See details.
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   evaluating and comparing multiple prediction vectors/models at once.
#' }
#' @param data A data.frame containing the observed outcome, assigned treatment,
#'   and necessary adjustment variables (confounders) for the evaluation of
#'   object.
#' @param outcome The outcome of interest within data. This could either be the
#'   name of a single numeric/logical column in data, or a Surv object for
#'   time-to-event data, e.g. Surv(time, status), if time and status are columns
#'   in data.
#' @param treatment_formula A formula which indicates the treatment/intervention variable
#'   (left hand side) and the adjustment variables (right hand side) in the
#'   data. E.g. A ~ L. The left hand side can be either a binary treatment
#'   (coded as 0/1 numeric, logical or factor) or a treatment with more than two
#'   categories (coded as a factor). The right hand side variables  are used to
#'   estimate the inverse probability of treatment weights (IPTW) with logistic/multinomial regression. The IPTW can
#'   also be specified directly as a vector using the iptw argument, in which case the right
#'   hand side of treatment_formula is ignored (the left hand side must still
#'   indicate the treatment, i.e. A ~ 1).
#' @param treatment_of_interest A treatment level under which the predictions
#'   should be evaluated.
#' @param metrics A character vector specifying which performance metrics to
#'   compute. Options are c("auc", "brier", “scaled_brier”, "oeratio",
#'   "calplot"). See details.
#' @param time_horizon For time to event data, the prediction horizon of
#'   interest.
#' @param cens_model Model for estimating inverse probability of censored
#'   weights (IPCW). Methods currently implemented are Kaplan-Meier ("KM") or
#'   Cox ("cox"), with censoring times derived from the Surv object specified
#'   under outcome, reversing the event indicator, see details. KM is only
#'   supported when the right hand side of cens_formula is 1.
#' @param cens_formula Model formula from which the right hand side is used in
#'   estimating the censoring probabilities. The left hand side must be left blank. E.g. ~ x1 + x2.
#' @param null_model If TRUE fits a model without covariates that estimates the same probability for all subjects in data. The model is
#'   fitted using the reweighted data in which all subjects 'counterfactually'
#'   received the treatment level of interest (using the IPTW, as estimated
#'   using the treatment_formula or as given by the iptw argument). For
#'   time-to-event outcomes, the null model is also fitted using the IPCW, as estimated using the cens_formula, or as given by the ipcw argument. The null_model can be used as reference (baseline) model.
#' @param bootstrap If this is an integer greater than 0, this indicates the
#'   number of bootstrap iterations, used to compute 95\% confidence intervals
#'   around the performance metrics based on percentiles of the bootstrap results.
#' @param bootstrap_progress if set to TRUE, print a progress bar indicating the
#'   progress of the bootstrap procedure.
#' @param iptw A numeric vector, containing the inverse probability of treatment
#'   weights. If iptw is not specified, these weights are  computed using the
#'   treatment_formula, but they can be specified directly via this argument. A
#'   user-defined function can also be specified, which takes as input 'data'
#'   and returns a numeric vector of IPTW weights. See details.
#' @param ipcw A numeric vector, containing the inverse probability of censoring
#'   weights at the time_horizon, or at a subject's event time, whichever happens
#'   first. For subjects who are censored before the time_horizon, the ipcw can be left at NA. If ipcw is not specified, these weights are computed using the
#'   cens_formula, but they can be specified directly via this argument. A
#'   user-defined function can also be specified, which takes as input 'data'
#'   and returns a numeric vector of IPCW weights. See details.
#' @param quiet If set to TRUE, don't print assumptions.
#' @param strip_ipt_models If set to TRUE (default), unnecessary components from
#'   the IPT- and IPC-model objects are not stored to save memory. Set to FALSE
#'   if you want to store the full IPT/IPC model objects.
#'
#' @returns An object of class `ip_score`, for which the `print()` and `plot()`
#' methods are implemented. The object is a nested list containing: \itemize{
#'   \item `$score`, contains the estimated predictive performance metrics.
#'   \item `$bootstrap`, if requested, the 95\% confidence intervals of the
#'   performance metrics, and the performance metrics for each individual
#'   bootstrap iteration.
#'   \item `$outcome`, the observed outcomes in data.
#'   \item `$treatment`, the observed treatment levels in data.
#'   \item `$predictions`, the predictions to be evaluated, i.e. the estimated probability
#'   of event under the intervention of interest for each subject.
#'   \item `$ipt`, method, model and inverse probability of treatment weights
#'   (IPTW). The IPTW are NA for subjects who did not receive the treatment level of interest.
#'   \item `$ipc`, method, model and inverse probability of censoring weights
#'   (IPCW). The IPCW are NA for subjects who were censored before the time_horizon.
#'   \item `$pseudopop`, binary vector indicating which subjects in data were re-weighted to create the pseudo-population. These are the subjects who were observed to receive the
#'   treatment level of interest and were not censored before the time_horizon, if applicable.
#'   }
#'   The print method summarizes the results and (if quiet = FALSE), prints the
#'   assumptions required for valid inference.
#'
#' @export
#'
#' @references Keogh RH, Van Geloven N. Prediction Under Interventions:
#'   Evaluation of Counterfactual Performance Using Longitudinal Observational
#'   Data. Epidemiology. 2024;35(3):329-339.
#'
#'   Boyer CB, Dahabreh IJ, Steingrimsson JA. Estimating and Evaluating
#'   Counterfactual Prediction Models. Statistics in Medicine.
#'   2025;44(23-24):e70287.
#'
#'   Pajouheshnia R, Peelen LM, Moons KGM, Reitsma JB, Groenwold RHH. Accounting
#'   for treatment use when validating a prognostic model: a simulation study.
#'   BMC Medical Research Methodology. 2017;17(1):103.
#'
#' @examples
#' n <- 1000
#'
#' data <- data.frame(L = rnorm(n), P = rnorm(n))
#' data$A <- rbinom(n, 1, plogis(data$L))
#' data$Y <- rbinom(n, 1, plogis(0.1 + 0.5*data$L + 0.7*data$P - 2*data$A))
#'
#' random <- runif(n, 0, 1)
#' model <- glm(Y ~ A + P, data = data, family = "binomial")
#'
#' score <- ip_score(
#'   object = list(random, model),
#'   data = data,
#'   outcome = Y,
#'   treatment_formula = A ~ L,
#'   treatment_of_interest = 0,
#' )
#' print(score)
#' plot(score)

ip_score <- function(object, data, outcome, treatment_formula,
                     treatment_of_interest,
                     metrics = c("auc", "brier", "scaled_brier", "oeratio", "calplot"),
                     time_horizon, cens_model = "KM", cens_formula = ~ 1,
                     null_model = TRUE,
                     bootstrap = 0, bootstrap_progress = TRUE,
                     iptw, ipcw, quiet = FALSE, strip_ipt_models = TRUE) {

  # checking inputs ---------------------------------------------------------


  check_missing(object)
  check_missing(data)
  check_missing(outcome)
  check_missing(treatment_formula)
  check_missing(treatment_of_interest)

  if (cens_model == "KM")
    stopifnot("censoring model must be ~ 1 if modeled via KM" =
                rhs_is_one(cens_formula))


  cens_model <- match.arg(cens_model, choices = c("cox", "KM"))

  # assert rhs(outcome_formula != 1) iff surv model AND!missing(iptw_weights)
  # assert longest surv time is longer than time horizon, to avoid annoying
  # weights

  object <- make_named_list(object, substitute(object))

  # do not allow bootstrap if iptw are given as fixed vector
  is_bootstrap_allowed(bootstrap, iptw, ipcw)

  # start gathering information required for the computation of metrics
  # in weighted pseudopop
  score_outcome <- extract_outcome(data, substitute(outcome), time_horizon)

  score_treatment <- extract_treatment(data, treatment_formula,
                                       treatment_of_interest)

  score_predictions <- get_predictions(object, data,
                                       score_treatment$treatment_column,
                                       score_treatment$treatment_of_interest,
                                       score_outcome$time_horizon)

  score_pseudopop <- get_pseudopop(score_outcome, score_treatment)

  score_ipt <- get_iptw(
    data = data,
    score_treatment = score_treatment,
    iptw = iptw,
    strip_model = strip_ipt_models
  )

  if (score_outcome$type == "survival") {
    cens_formula <- combine_censoring_formula(cens_formula, substitute(outcome))
    score_ipc <- get_ipcw(cens_formula, data, cens_model, time_horizon, ipcw)
  } else {
    score_ipc <- NULL
  }

  if (null_model) {
    score_predictions <- fit_null(score_pseudopop, score_outcome,
                                  score_predictions, score_ipt, score_ipc)
  }

  # make object
  ip_object <- construct_ip_object(score_outcome, score_treatment,
                                   score_predictions, score_pseudopop,
                                   score_ipt, score_ipc,
                                   metrics)

  # compute metrics
  ip_object <- compute_metrics(ip_object)

  # do bootstrap
  if (bootstrap > 0) {
    matchcall <- match.call()
    call_env <- parent.frame()

    bs <- bootstrap(ip_object, matchcall, call_env, bootstrap, bootstrap_progress)
    ip_object <- add_to_ip_object(ip_object, "bootstrap", bs, after = 1)
    ip_object <- add_to_ip_object(ip_object, "bootstrap_iterations", bootstrap)
  }

  ip_object <- add_to_ip_object(ip_object, "quiet", quiet)
  return(ip_object)
}

handle_specified_ip <- function(iptcw, data, type = "iptw") {
  # if the user specified something in the iptw/ipcw argument, e.g. a vector of
  # weights or a function:
  if (is.vector(iptcw, mode = "numeric")) {

    method <- "weights manually specified"
    ipw <- iptcw
  } else if (is.function(iptcw)) {

    method <- "weights specified via function"

    ipw <- do.call(iptcw, args = list("data" = data))

    if ( ! (is.vector(ipw, mode = "numeric") && length(ipw) == nrow(data))) {

      stop("function specified in ", type, " did not return a numeric vector ",
           "of length data")
    }
  } else {

    stop("argument ", type, " must be missing, a numeric vector of weights, ",
         "or a function.")

  }

  return(list("method" = method, "ipw" = ipw))

}

get_iptw <- function(data, score_treatment, iptw,
                     only_weights = FALSE, strip_model = TRUE) {
  ipt <- list()
  # if user specified weights themselves:
  if (!missing(iptw)) {
    manualiptw <- handle_specified_ip(iptw, data)
    ipt$method <- manualiptw$method
    iptw <- manualiptw$ipw
  } else {

    # else we compute the weights ourselves.
    trt_form <- score_treatment$propensity_formula
    ipt$confounders <- all.vars(trt_form)[-1]
    ipt$propensity_formula <- trt_form
    iptw_object <- ipt_weights(
      data = data,
      propensity_formula = trt_form,
      treatment_of_interest = score_treatment$treatment_of_interest,
      type = score_treatment$type,
      strip_model = strip_model
    )
    ipt$model <- iptw_object$model
    ipt$method <- iptw_object$method
    iptw <- iptw_object$weights

  }
  ipt$weights <- as.vector(iptw)
  if (only_weights == TRUE) {
    return(list("weights" = ipt$weights))
  } else {
    return(ipt)
  }
}

combine_censoring_formula <- function(cens_formula, outcome) {
  stats::update.formula(
    old = cens_formula,
    new = stats::as.formula(call("~", outcome, quote(.)))
  )
}

compute_metrics <- function(ip_object) {
  metrics <- list()

  if (ip_object$outcome$type == "survival") {
    weights <- ip_object$ipt$weights * ip_object$ipc$weights
    outcome <- ip_object$outcome$status_at_horizon
  } else {
    weights <- ip_object$ipt$weights
    outcome <- ip_object$outcome$observed
  }

  for (m in ip_object$metrics) {
    metrics[[m]] <- sapply(
      X = ip_object$predictions,
      FUN = function(x) {
        cf_metric(
          m,
          obs_outcome = outcome,
          cf_pred = x,
          pseudo_i = ip_object$pseudopop$ids,
          ipw = weights
        )
      }
    )
  }

  ip_object <- add_to_ip_object(ip_object, "score", metrics, after = 0)

  return(ip_object)
}

get_pseudopop <- function(score_outcome, score_treatment) {
  pseudopop_list <- list()

  correct_trt <- with(score_treatment, observed == treatment_of_interest)

  if (score_outcome$type == "survival") {
    uncensored <- with(
      score_outcome,
      observed[, 1] >= time_horizon | observed[, 2] == 1
    )
    pseudopop <- correct_trt & uncensored

    pseudopop_list$ids <- pseudopop
    pseudopop_list$correct_trt <- correct_trt
    pseudopop_list$uncensored <- uncensored

  } else {
    pseudopop_list$ids <- correct_trt
    pseudopop_list$correct_trt <- correct_trt
  }
  return(pseudopop_list)
}

construct_ip_object <- function(outcome, treatment, predictions, pseudopop,
                                ipt, ipc, metrics) {

  ip_object <- list(
    "outcome" = outcome,
    "treatment" = treatment,
    "predictions" = predictions,
    "pseudopop" = pseudopop,
    "metrics" = metrics,
    "ipt" = ipt
  )
  if (ip_object$outcome$type == "survival") {
    ip_object$ipc = ipc
  }

  class(ip_object) <- "ip_score"
  return(ip_object)
}

add_to_ip_object <- function(ip_object, name, value, after = length(ip_object)) {
  old_class <- class(ip_object)

  ip_object <- append(ip_object, list(value), after = after)
  names(ip_object)[after+1] <- name
  class(ip_object) <- old_class

  return(ip_object)
}

extract_outcome <- function(data, outcome, time_horizon) {
  # attempt to extract the outcome from the data, and perform various sanity
  # checks
  y <- tryCatch(
    eval(outcome, envir = as.list(data), enclos = parent.frame()),
    error = function(e) {
      stop(sprintf("Outcome %s not found in data", deparse(outcome)),
           call. = FALSE)
    }
  )

  if (!( ((is.numeric(y) || is.logical(y)) && is.vector(y)) ||
         inherits(y, "Surv") )) {
    stop("Outcome must be a numeric vector or a Surv object", call. = FALSE)
  }

  if (length(y) != nrow(data)) {
    stop("Outcome must be of length nrow(data)", call. = FALSE)
  }

  if (inherits(y, "Surv")) {
    time  <- y[, 1]
    event <- y[, 2]

    if (any(!is.finite(time)) || any(time < 0)) {
      stop("Survival times must be finite & nonnegative", call. = FALSE)
    }

    if (!all(event %in% c(0, 1))) {
      stop("Event indicator must be binary (0/1)", call. = FALSE)
    }

  } else {
    if (!all(y %in% c(0, 1))) {
      stop("Outcome must be binary (0/1)", call. = FALSE)
    }
  }

  outcome_list <- list("observed" = y)

  if (inherits(outcome_list$observed, "Surv")) {
    outcome_list$type <- "survival"
    outcome_list$time_horizon <- time_horizon
    outcome_list$status_at_horizon <- ifelse(
      test = outcome_list$observed[, 1] <= time_horizon,
      yes = outcome_list$observed[, 2],
      no = FALSE
    )
  } else {
    outcome_list$type <- "binary"
  }
  outcome_list
}

extract_treatment <- function(data, treatment_formula, treatment_of_interest) {
  trt_list <- list()
  trt_list$treatment_column <- treatment_formula[[2]]
  trt_list$observed <- extract_lhs(data, treatment_formula)
  trt_list$treatment_of_interest <- treatment_of_interest
  trt_list$propensity_formula <- treatment_formula

  n_trt <- length(unique(trt_list$observed))

  if (n_trt == 2) {
    trt_list$type <- "binary"
  } else if (n_trt >= 3) {
    stopifnot(
      "More than 2 treatment options found in data, but not a factor variable." =
        is.factor(trt_list$observed)
    )
    trt_list$type <- "categorical"
  } else {
    if (trt_list$treatment_column != "ipscore_fake_trt") {
      stop("Only 1 treatment option found in data. Must have at least 2 options.")
    }
  }
  if (!is.na(treatment_of_interest)) {
    stopifnot("Specified treatment_of_interest value does not appear in data" =
                treatment_of_interest %in% trt_list$observed)
  }

  trt_list
}

get_predictions <- function(object, data, treatment_column,
                            treatment_of_interest, time_horizon) {

  predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && (is.null(dim(x)) || length(dim(x)) == 1)) {
        x <- as.vector(x)
        stopifnot("Predictions must be of length nrow(data)" =
                    length(x) == nrow(data))
        stopifnot("Predictions must be in interval [0,1]" =
                    all(x >= 0) && all(x <= 1))
        x # user supplied risk predictions
      } else {
        predict_CF(
          x,
          data,
          treatment_column,
          treatment_of_interest,
          time_horizon
        )
      }
    }
  )
  predictions
}

get_ipcw <- function(cens_formula, data, cens_model, time_horizon,
                     ipcw, only_weights = FALSE, strip_ipt_models = TRUE) {
  ipc <- list()

  # if user specified weights themselves:
  if (!missing(ipcw)) {
    manualiptw <- handle_specified_ip(ipcw, data, type = "ipcw")
    ipc$method <- manualiptw$method
    ipcw <- manualiptw$ipw
  } else {

    ipc$method <- cens_model
    ipc$cens_formula <- cens_formula
    ipc_object <- ipc_weights(data, ipc$cens_formula,
                              cens_model, time_horizon, strip_ipt_models)
    ipcw <- ipc_object$weights
    ipc$model <- ipc_object$model

    attr(ipc$cens_formula, ".Environment") <- NULL

  }
  ipc$weights <- ipcw

  if (only_weights == TRUE) {
    return(list("weights" = ipc$weights))
  } else {
    return(ipc)
  }
}

fit_null <- function(score_pseudopop, score_outcome, score_predictions,
                     score_ipt, score_ipc) {
  # fit a null on the pseudo-population that received treatment of
  # interest and remained uncensored. Add it to score_predictions.
  pseudo_ids <- score_pseudopop[[1]]
  n <- length(score_outcome$observed)
  if (score_outcome$type == "binary") {
    outcomes <- score_outcome$observed[pseudo_ids]
    weights <- score_ipt$weights[pseudo_ids]
  } else {
    outcomes <- score_outcome$status_at_horizon[pseudo_ids]
    weights <- score_ipt$weights[pseudo_ids]*score_ipc$weights[pseudo_ids]
  }

  null_prediction <- stats::weighted.mean(outcomes, weights)
  null_preds <- rep(null_prediction, n)

  score_predictions <- c(
    list("null model" = null_preds),
    score_predictions
  )
  return(score_predictions)
}
