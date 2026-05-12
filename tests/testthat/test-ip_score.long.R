test_that("ipscore long works", {
  n_dev <- 1000
  n_val <- 200000

  df_dev <- generate_long_data_cox(n_dev, seed = 1)
  df_dev_long <- make_dev_long(df_dev)
  iptw <- ipt_weights(df_dev_long, A ~ L * A_lag_1)$weights

  coxmsm <- fit_long_cox_model(data_long = df_dev_long, iptw)
  badmodel <- glm(status ~ A0 + L0, data = df_dev, family = "binomial")

  df_val <- generate_long_data_cox(n_val, seed = 3)
  df_cf0 <- generate_long_data_cox(n_val, seed = 3, Ai = function(i) rep(0, n_val))
  df_cf1 <- generate_long_data_cox(n_val, seed = 3, Ai = function(i) rep(1, n_val))

  df_val_outcome <- df_val[, c("id", "time", "status")]
  df_val_long <- wide_to_long(df_val, "id", list(A = paste0("A", 0:4),
                                                 L = paste0("L", 0:4)),
                              0:4, df_val$time)
  df_val_long <- add_lag_terms(df_val_long, "A")


  risk_0 <- risk_under_0(coxmsm, 5, df_val$L0)
  risk_1 <- risk_under_1(coxmsm, 5, df_val$L0)
  risk_bad_0 <- predict_CF(badmodel, df_val, "A0", 0)
  risk_bad_1 <- predict_CF(badmodel, df_val, "A0", 1)
  risk_random <- runif(n_val)
  always_0 <- rep(0, n_val)
  always_1 <- rep(1, n_val)
  always_truth_under_0 <- df_cf0$status
  always_wrong_under_0 <- 1-df_cf0$status
  always_truth_under_1 <- df_cf1$status
  always_wrong_under_1 <- 1-df_cf1$status


  metrics <- c("auc", "brier", "oeratio", "scaled_brier")

  models <- list(risk_0, risk_1, risk_bad_0, risk_bad_1, risk_random,
                 always_0, always_1, always_truth_under_0, always_wrong_under_0,
                 always_truth_under_1, always_wrong_under_1)

  score0 <- ip_score_long(
    probabilities = models,
    data_outcome = df_val_outcome,
    data_long = df_val_long,
    time_horizon = 5,
    treatment_formula = A ~ (A_lag_1 * L),
    treatment_of_interest = rep(0, 5),
    null_model = TRUE,
    metrics = metrics
  )
  score0_true <- observed_score(models, df_cf0, df_cf0$status, metrics, TRUE)

  score1 <- ip_score_long(
    probabilities = models,
    data_outcome = df_val_outcome,
    data_long = df_val_long,
    time_horizon = 5,
    treatment_formula = A ~ (A_lag_1 * L),
    treatment_of_interest = rep(1, 5),
    null_model = TRUE,
    metrics = metrics
  )
  score1_true <- observed_score(models, df_cf1, df_cf1$status, metrics, TRUE)

  expect_equal(score0$score, score0_true$score, tolerance = 0.01)
  expect_equal(score1$score, score1_true$score, tolerance = 0.01)

})
