ipt_weights <- function(data, propensity_formula, treatment_of_interest,
                        strip_model = TRUE) {

  trt_variable <- all.vars(propensity_formula)[1]
  propensity_model <- stats::glm(propensity_formula, family = "binomial", data,
                                 x = FALSE, y = FALSE, model = FALSE)
  prop_score <- unname(stats::predict(propensity_model, type = "response"))

  if (strip_model) {
    propensity_model <- strip_glm(propensity_model)
  }

  prob_trt <- ifelse(data[[trt_variable]] == 1, prop_score, 1 - prop_score)

  weights <- 1 / prob_trt

  if (!missing(treatment_of_interest)) {
    weights[data[[trt_variable]] != treatment_of_interest] <- NA
  }

  list(
    model = propensity_model,
    weights = weights
  )
}


