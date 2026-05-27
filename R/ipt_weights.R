ipt_weights <- function(data, propensity_formula, treatment_of_interest = NA,
                        type = "binary", strip_model = TRUE) {

  trt_variable <- all.vars(propensity_formula)[1]

  model_prob <- switch(
    type,
    "binary" = ipt_model_prob_binary(data, propensity_formula, trt_variable),
    "categorical" = ipt_model_prob_categorical(data, propensity_formula, trt_variable),
    stop(paste0("iptw type ", type, " not implemented!"))
  )

  if (strip_model) {
    model_prob[[1]] <- strip_glm(model_prob[[1]])
  }

  weights <- 1 / model_prob[[2]]

  # set weights of patients outside pseudopopulation to NA
  if (!is.na(treatment_of_interest)) {
    weights[data[[trt_variable]] != treatment_of_interest] <- NA
  }

  list(
    model = model_prob[[1]],
    weights = as.vector(weights),
    method = model_prob[[3]]
  )
}


ipt_model_prob_binary <- function(data, propensity_formula, trt_variable) {
  propensity_model <- stats::glm(propensity_formula, family = "binomial", data,
                                 x = FALSE, y = FALSE, model = FALSE)
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

  list(propensity_model, prob_trt, "binomial glm")
}

ipt_model_prob_categorical <- function(data, propensity_formula, trt_variable) {

  propensity_model <- nnet::multinom(propensity_formula, data = data)
  probability_of_treatments <- stats::predict(propensity_model, data,
                                              type = "probs")

  # select probability corresponding to observed treatment
  prob_trt <- probability_of_treatments[
    cbind(
      1:nrow(probability_of_treatments),
      data[[trt_variable]]
    )
  ]
  list(propensity_model, prob_trt, "multinomial")
}
