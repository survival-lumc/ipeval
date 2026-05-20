ipt_weights <- function(data, propensity_formula, treatment_of_interest,
                        type = "binomial") {

  trt_variable <- all.vars(propensity_formula)[1]

  if (type == "binomial") {
    propensity_model <- stats::glm(propensity_formula, family = "binomial", data)
    prop_score <- unname(stats::predict(propensity_model, type = "response"))

    if (is.factor(data[[trt_variable]])) {
      # if trt variable is a factor, predict(propensity_model) estimates
      # probability of the second label. if binary, predict(propensity_model)
      # estimates probability of outcome 1.
      predicted_level <- levels(data[[trt_variable]])[[2]]
    } else {
      predicted_level <- 1
    }
    # probability of getting the treatment that the patient was assigned
    prob_trt <- ifelse(
      data[[trt_variable]] == predicted_level,
      prop_score,
      1 - prop_score
    )
  } else if (type == "categorical") {
    propensity_model <- nnet::multinom(propensity_formula, data = data)
    probability_of_treatments <- stats::predict(propensity_model, data,
                                                type = "probs")
    prob_trt <- probability_of_treatments[, data[[trt_variable]]]
  }

  weights <- 1 / prob_trt

  # set weights of patients outside pseudopopulation to NA
  if (!missing(treatment_of_interest)) {
    weights[data[[trt_variable]] != treatment_of_interest] <- NA
  }

  list(
    model = propensity_model,
    weights = weights
  )
}

