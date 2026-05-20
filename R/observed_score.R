#' Performance in observed dataset
#'
#' This function computes the performance of the predictions in the given data,
#' which may contain a mix of treated and untreated subjects. It exists only to
#' demonstrate the difference between 'normal' performance and counterfactual
#' performance. It is not user friendly and should not be relied on. It does not
#' support time-to-event data.
#'
#' @param object One of the following three options to be validated:
#' \itemize{
#'   \item a numeric vector, corresponding to risk predictions
#'   \item a glm model
#'   \item a (named) list, with one or more of the previous 2 options, for
#'   validating and comparing multiple models at once.
#' }
#' @param data A data.frame containing the observed outcome.
#' @param outcome The outcome, to be evaluated within data. This should be
#'   the name of a numeric column in data.
#' @param metrics A character vector specifying which performance metrics to be
#'   computed. Options are c("auc", "brier", "oeratio", "calplot").
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
                    metrics = c("auc", "brier", "oeratio", "calplot")) {

  # make a list of risk predictions
  object <- make_list_if_not_list(object)
  score_predictions <- lapply(
    X = object,
    FUN = function(x) {
      if (is.numeric(x) && is.null(dim(x))) {
        x # user supplied risk predictions
      } else {
        stats::predict(x, newdata = data, type = "response")
      }
    }
  )

  # performance of observed data, can be computed by setting everyone's treatment
  # to the treatment of interest. This way all subjects are used
  score_trt <- list(
    "observed" = rep(1, nrow(data)),
    "treatment_of_interest" = 1
  )

  ip_object <- construct_ip_object(
    outcome = extract_outcome(data, substitute(outcome)),
    treatment = score_trt,
    predictions = score_predictions,
    pseudopop = list(ids = rep(TRUE, nrow(data))),
    ipt = list("weights" = rep(1, nrow(data))),
    ipc = NULL,
    metrics = metrics
  )
  ip_object <- add_to_ip_object(ip_object, "quiet", TRUE)
  ip_object <- compute_metrics(ip_object)

  return(ip_object)
}
