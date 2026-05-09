
n <- 10000

run_test <- function(seed) {
  set.seed(seed)
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  model <- suppressWarnings(
    glm(
      Y ~ A + P,
      family = "binomial",
      data = data,
      weights = ipt_weights(data, A ~ L)$weights
    )
  )

  Y0_predicted <- predict_CF(model, data, "A", 0)

  ip_score <- ip_score(
    data = data,
    object = Y0_predicted,
    outcome = Y,
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    null_model = TRUE,
    metrics = c("brier", "scaled_brier")
  )


  score <- riskRegression::Score(
    list(Y0_predicted),
    formula = Y0 ~ 1,
    data = data,
    null.model = T,
    metrics = c("brier")
  )

  # score2 <- observed_score(
  #   Y0_predicted, data, Y0, metrics = c("brier", "scaled_brier")
  # )

  score_brier_model <- score$Brier$score[[2, "Brier"]]
  score_brier_null <- score$Brier$score[[1, "Brier"]]
  scaled_brier <- (1 - score_brier_model/score_brier_null)*100

  # expect_equal(score2$score$brier[[1]], score_brier_model)
  # expect_equal(score2$score$scaled_brier[[1]], scaled_brier)
  #
  return(c(
    "null" = ip_score$predictions[[1]][[1]],
    "brier_null" = ip_score$score$brier[[1]],
    "brier_model" = ip_score$score$brier[[2]],
    "scaled_brier_model" = ip_score$score$scaled_brier[[2]],
    "true_null" = mean(data$Y0),
    "true_brier_null" = score_brier_null,
    "true_brier_model" = score_brier_model,
    "true_scaled_brier_model" = scaled_brier
  ))
}

results <- lapply_progress(as.list(1:500), run_test, "")
results <- as.data.frame(do.call(rbind, results))

results$bias_null_model <- results$null - results$true_null
results$bias_brier_null <- results$brier_null - results$true_brier_null
results$bias_brier_model <- results$brier_model - results$true_brier_model
results$bias_scaled <- results$scaled_brier_model - results$true_scaled_brier_model

se <- function(x) sd(x) / sqrt(length(x))

apply(results, 2, mean)
apply(results, 2, se)

t.test(results$bias_null_model)
t.test(results$bias_brier_null) # huh?
t.test(results$bias_brier_model)
t.test(results$bias_scaled)
