#' Performance in observed dataset
#'
#' This function computes the performance of the predictions in the given data,
#' which may contain a mix of treated and untreated subjects. It exists only to
#' demonstrate the difference between 'normal' performance and counterfactual
#' performance. It is not user friendly and should not be relied on. Consider
#' using riskRegression::Score() as an alternative.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions
#'   \item a glm model
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome.
#' @param outcome The outcome, to be evaluated within data. This could either be
#'   the name of a numeric/logical column in data, or a Surv object for
#'   time-to-event data, e.g. Surv(time, status), if time and status are columns
#'   in data.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot").
#' @param time_horizon For time to event data, the prediction horizon of
#'   interest.
#' @param cens_model Model for estimating inverse probability of censored
#'   weights (IPCW). Methods currently implemented are Kaplan-Meier ("KM") or
#'   Cox ("cox"), both applied to the censored times. KM is only supported when
#'   the right hand side of cens_formula is 1.
#' @param cens_formula Formula for which the r.h.s. determines the censoring
#'   probabilities. I.e. ~ x1 + x2.
#' @param null_model If TRUE fit a risk prediction model which ignores the
#'   covariates and predicts the same value for all subjects. For time-to-event
#'   outcomes, the subjects are 'counterfactually' uncensored (using the
#'   IPCW, as estimated using the cens_formula, or as given by the ipcw
#'   argument).
#' @param ipcw A numeric vector, containing the inverse probability of censor
#'   weights. These are normally computed using the cens_formula, but they can
#'   be specified directly via this argument.
#'
#' @returns Performance metrics in the observed dataset.
#' @export
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
#' observed_score(
#'   object = list("ran" = random, "mod" = model),
#'   data = data,
#'   outcome = Y,
#'   metrics = c("auc", "brier", "oeratio")
#' )

observed_score <- function(object, data, outcome,
  metrics = c("auc", "brier", "oeratio", "calplot"),
  time_horizon, cens_model = "KM", cens_formula = ~ 1,
  null_model = TRUE, ipcw) {

  # observed score can be achieved by calling ip_score with some fake treatment
  # value that everybody has, and with iptw = 1 for everybody

  if ("ipscore_fake_trt" %in% names(data))
    stop("column name ipscore_fake_trt is reserved for internal use.")

  data$ipscore_fake_trt <- 1
  data2 <- data
  treatment_of_interest <- 1
  treatment_formula <- ipscore_fake_trt ~ 1

  cl <- match.call()
  cl[[1]] <- quote(ip_score)

  cl$data <- data2
  cl$treatment_formula <- quote(ipscore_fake_trt ~ 1)
  cl$treatment_of_interest <- 1
  cl$stable_iptw <- FALSE
  cl$iptw <- rep(1, nrow(data))

  observed_score <- eval.parent(cl)

  # remove artifacts from ip_score that are not important for observed score
  observed_score$treatment <- NULL
  observed_score$pseudopop <- NULL
  observed_score$ipt <- NULL
  observed_score$quiet <- NULL
  class(observed_score) <- NULL

  return(observed_score)
}
