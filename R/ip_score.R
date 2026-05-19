#' Counterfactual validation score
#'
#' Estimates the predictive performance of predictions under interventions, by
#' forming a weighted pseudopopulation in which every subject was assigned the
#' treatment of interest.
#'
#' When supplying a glm or coxph model, the model should be able to estimate
#' risks under the intervention of interest. This could be done in two ways: the
#' model does not have a treatment covariate, and always estimates the risk
#' under intervention of interest. Alternatively, the model has a covariate for
#' treatment. This function then automatically estimates the risk under the
#' treatment of interest for all subjects, even if they were assigned
#' alternative treatment.
#'
#' All performance metrics are computed on the weighted population in which
#' every subject was counterfactually assigned treatment of interest. \code{auc}
#' is area under the (ROC) curve. Brier score is defined as \code{1 / sum(iptw)
#' sum(predictions_i - outcome_i)^2}. Scaled brier score is also available
#' (\code{metrics = "scaled_brier"}). For the O/E ratio, the numerator
#' (observed) is the weighted fraction of events in the pseudopopulation, and
#' the denominator (expected) is the unweighted mean of risk estimates of the
#' original unweighted population. The \code{calplot} option generates a
#' calibration plot, with default 8 knots. More/less knots can be specified by
#' appending calplot with a number indicating the number of knots, e.g.
#' \code{metrics = calplot10} for 10 knots.
#'
#' Stabilized IPT-weigths are computed by estimating a null model for treatment.
#' E.g. weights are \code{P(A = a) / P(A = a | L = l)}, if the given
#' treatment_formula is \code{A ~ L}.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions under
#'   intervention of interest.
#'   \item a glm or coxph model, capable of estimating risks under intervention
#'   of interest. See details.
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome, assigned treatment,
#'   and necessary confounders for the validation of object.
#' @param outcome The outcome, to be evaluated within data. This could either be
#'   the name of a numeric/logical column in data, or a Surv object for
#'   time-to-event data, e.g. Surv(time, status), if time and status are columns
#'   in data.
#' @param treatment_formula A formula which identifies the treatment (left hand
#'   side) and the confounders (right hand side) in the data. E.g. A ~ L. The
#'   confounders are used to estimate the inverse probability of treatment
#'   weights (IPTW) model. The IPTW can also be specified themselves using the
#'   iptw argument, in which case the right hand side of this formula is
#'   ignored.
#' @param treatment_of_interest A treatment level for which the counterfactual
#'   perormance measures should be evaluated.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot"). See details.
#' @param time_horizon For time to event data, the prediction horizon of
#'   interest.
#' @param cens_model Model for estimating inverse probability of censored
#'   weights (IPCW). Methods currently implemented are Kaplan-Meier ("KM") or
#'   Cox ("cox"), both applied to the censored times. KM is only supported when
#'   the right hand side of cens_formula is 1.
#' @param cens_formula Formula for which the r.h.s. determines the censoring
#'   probabilities. I.e. ~ x1 + x2.
#' @param null_model If TRUE fit a risk prediction model which ignores the
#'   covariates and predicts the same value for all subjects. The model is
#'   fitted using the data in which all subjects are counterfactually assigned
#'   the treatment of interest (using the IPTW, as estimated using the
#'   treatment_formula or as given by the iptw argument). For time-to-event
#'   outcomes, the subjects are also counterfactually uncensored (using the
#'   IPCW, as estimated using the cens_formula, or as given by the ipcw
#'   argument).
#' @param stable_iptw if TRUE, estimate stabilized IPT-weights. See details.
#' @param bootstrap If this is an integer greater than 0, this indicates the
#'   number of bootstrap iterations, to compute 95\% confidence intervals around
#'   the performance metrics.
#' @param bootstrap_progress if set to TRUE, print a progress bar indicating the
#'   progress of bootstrap procedure.
#' @param iptw A numeric vector, containing the inverse probability of treatment
#'   weights. These are normally computed using the treatment_formula, but they
#'   can be specified directly via this argument. If specified via this
#'   argument, bootstrap is not possible.
#' @param ipcw A numeric vector, containing the inverse probability of censor
#'   weights. These are normally computed using the cens_formula, but they can
#'   be specified directly via this argument. If specified via this argument,
#'   bootstrap is not possible.
#' @param quiet If set to TRUE, don't print assumptions.
#'
#' @returns A list with performance metrics, of class 'ip_score", for which the
#' print and plot methods are implemented.
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
#' naive_perfect <- data$Y
#'
#' score <- ip_score(
#'   object = list("ran" = random, "mod" = model, "per" = naive_perfect),
#'   data = data,
#'   outcome = Y,
#'   treatment_formula = A ~ L,
#'   treatment_of_interest = 0,
#' )
#' print(score)
#' plot(score)

ip_score <- function(object, data, outcome, treatment_formula,
                     treatment_of_interest,
                     metrics = c("auc", "brier", "oeratio", "calplot"),
                     time_horizon, cens_model = "KM", cens_formula = ~ 1,
                     null_model = TRUE, stable_iptw = FALSE,
                     bootstrap = 0, bootstrap_progress = TRUE,
                     iptw, ipcw, quiet = FALSE) {


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

  if (bootstrap != 0)
    stopifnot("can't bootstrap if iptw are given" = missing(iptw))

  score_outcome <- extract_outcome(data, substitute(outcome), time_horizon)

  score_treatment <- extract_treatment(data, treatment_formula,
                                       treatment_of_interest)

  score_predictions <- get_predictions(object, data,
                                       score_treatment$treatment_column,
                                       score_treatment$treatment_of_interest,
                                       score_outcome$time_horizon)

  score_ipt <- get_iptw(treatment_formula, data, stable_iptw, iptw)

  if (score_outcome$type == "survival") {
    cens_formula <- combine_censoring_formula(cens_formula, substitute(outcome))
    score_ipc <- get_ipcw(cens_formula, data, cens_model, time_horizon, ipcw)
  } else {
    score_ipc <- NULL
  }

  if (null_model) {
    score_predictions <- fit_null(score_treatment, score_outcome,
                                  score_predictions, score_ipt, score_ipc)
  }

  # make object
  ip_object <- construct_ip_object(score_outcome, score_treatment,
                                   score_predictions, score_ipt, score_ipc,
                                   metrics)

  # compute metrics
  ip_object <- compute_metrics(ip_object)

  # do bootstrap
  if (bootstrap > 0) {
    bs <- bootstrap(data, ip_object, bootstrap, bootstrap_progress)
    ip_object <- add_to_ip_object(ip_object, "bootstrap", bs, after = 1)
    ip_object <- add_to_ip_object(ip_object, "bootstrap_iterations", bootstrap)
  }

  ip_object <- add_to_ip_object(ip_object, "quiet", quiet)
  return(ip_object)
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
          obs_trt = ip_object$treatment$observed,
          cf_pred = x,
          cf_trt = ip_object$treatment$treatment_of_interest,
          ipw = weights
        )
      }
    )
  }

  ip_object <- add_to_ip_object(ip_object, "score", metrics, after = 0)

  return(ip_object)
}

# TODO: documentation page for structure of ip_object?
construct_ip_object <- function(outcome, treatment, predictions, ipt, ipc,
                                metrics) {

  ip_object <- list(
    "outcome" = outcome,
    "treatment" = treatment,
    "predictions" = predictions,
    "metrics" = metrics,
    "ipt" = ipt
  )
  if (ip_object$outcome$type == "survival") {
    ip_object$ipc = ipc
  }
  ip_object$correct_trt <- with(
    ip_object$treatment,
    observed == treatment_of_interest
  )
  if (ip_object$outcome$type == "survival") {
    ip_object$uncensored <- with(
      ip_object$outcome,
      ! ( observed[, 1] >= time_horizon | observed[, 2] == 1 )
    )
    ip_object$pseudopop <- ip_object$correct_trt & ip_object$uncensored
  } else {
    ip_object$pseudopop <- ip_object$correct_trt
  }

  class(ip_object) <- "ip_score"
  return(ip_object)
}

add_to_ip_object <- function(ip_object, name, value, after = length(ip_object)) {
  ip_object <- append(ip_object, list(value), after = after)
  names(ip_object)[after+1] <- name
  class(ip_object) <- "ip_score"
  return(ip_object)
}

name_unnamed_list <- function(x) {
  # give names, if not named
  sapply(
    1:length(x),
    function(i)

      if (is.null(names(x)[i]) || names(x)[i] == "") {
        paste0("model.", i)
      } else {
        names(x[i])
      }
  )
}

make_list_if_not_list <- function(x) {
  if (!("list" %in% class(x)))
    x <- list(x)
  names(x) <- name_unnamed_list(x)
  x
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
  stopifnot("Treatment is not binary" =
              setequal(unique(trt_list$observed), c(0,1)))
  stopifnot("Treatment_of_interest must be either 0 or 1" =
              treatment_of_interest == 0 || treatment_of_interest == 1)
  trt_list
}

get_predictions <- function(object, data, treatment_column,
                            treatment_of_interest, time_horizon) {
  # make a list of risk predictions
  object <- make_list_if_not_list(object)
  predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && is.null(dim(x))) {
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

get_iptw <- function(treatment_formula, data, stable_iptw, iptw,
                     only_weights = FALSE) {
  ipt <- list()
  ipt$method = "weights manually specified"

  if (missing(iptw)) {
    ipt$method <- "binomial glm"
    ipt$confounders <- all.vars(treatment_formula)[-1]
    ipt$propensity_formula <- treatment_formula
    iptw_object <- ipt_weights(data, treatment_formula)
    ipt$model <- iptw_object$model
    iptw <- iptw_object$weights

    if (stable_iptw == TRUE) {
      ipt$method <- "stabilized weights"
      stable_treatment_formula <-
        stats::update.formula(treatment_formula, . ~ 1)
      sipt_object <- ipt_weights(data, stable_treatment_formula)
      iptw <- 1/sipt_object$weights * iptw
      ipt$stable_model <- sipt_object$model
    }
  }
  ipt$weights <- iptw

  if (only_weights == TRUE) {
    return(list("weights" = ipt$weights))
  } else {
    return(ipt)
  }
}

get_ipcw <- function(cens_formula, data, cens_model, time_horizon,
                     ipcw, only_weights = FALSE) {
  ipc <- list()
  ipc$method <- "weights manually specified"
  if (missing(ipcw)) {
    ipc$method <- cens_model
    ipc$cens_formula <- cens_formula
    ipc_object <- ipc_weights(data, ipc$cens_formula,
                       cens_model, time_horizon)
    ipcw <- ipc_object$weights
    ipc$model <- ipc_object$model
  }
  ipc$weights <- ipcw

  if (only_weights == TRUE) {
    return(list("weights" = ipc$weights))
  } else {
    return(ipc)
  }
}

fit_null <- function(score_treatment, score_outcome, score_predictions,
                     score_ipt, score_ipc) {
  # fit a null on the pseudo-population that received treatment of
  # interest, and add it to the score_predictions
  pseudo_ids <- score_treatment$observed == score_treatment$treatment_of_interest
  n <- length(score_outcome$observed)
  if (score_outcome$type == "binary") {
    null_prediction <- stats::weighted.mean(
      score_outcome$observed[pseudo_ids],
      score_ipt$weights[pseudo_ids]
    )
  } else {
    # fit null model in a pseudopopuluation where everyone had treatment of
    # interest and remained uncensored
    uncensor_ids <- score_ipc$weights != 0
    cf_ids <- pseudo_ids & uncensor_ids
    null_prediction <- stats::weighted.mean(
      score_outcome$status_at_horizon[cf_ids],
      score_ipt$weights[cf_ids]*score_ipc$weights[cf_ids]
    )
  }
  null_preds <- rep(null_prediction, n)

  score_predictions <- c(
    list("null model" = null_preds),
    score_predictions
  )
  return(score_predictions)
}
