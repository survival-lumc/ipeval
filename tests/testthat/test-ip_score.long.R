test_that("ipscore long works", {
  n_dev <- 1000
  n_val <- 300000

  df_dev <- generate_long_data_cox(n_dev, seed = 1)
  df_dev_long <- make_dev_long(df_dev)
  iptw <- ipt_weights(df_dev_long, A ~ L * A_lag_1)$weights

  coxmsm <- fit_long_cox_model(data_long = df_dev_long, iptw)

  df_val <- generate_long_data_cox(n_val, seed = 2)
  df_cf0 <- generate_long_data_cox(n_val, seed = 2, Ai = function(i) 0)
  df_cf1 <- generate_long_data_cox(n_val, seed = 2, Ai = function(i) 1)

  df_val_outcome <- df_val[, c("id", "time", "status")]
  df_val_long <- wide_to_long(df_val, "id", list(A = paste0("A", 0:4),
                                                 L = paste0("L", 0:4)),
                              0:4, df_val$time)
  df_val_long <- add_lag_terms(df_val_long, "A")


  risk_untreated <- risk_under_0(coxmsm, 5, df_val$L0)
  risk_treated <- risk_under_1(coxmsm, 5, df_val$L0)

  metrics <- c("auc", "brier", "oeratio")

  score0 <- ip_score_long(
    probabilities = risk_untreated,
    data_outcome = df_val_outcome,
    data_long = df_val_long,
    time_horizon = 5,
    treatment_formula = A ~ (A_lag_1 * L),
    treatment_of_interest = rep(0, 5),
    null_model = FALSE,
    metrics = metrics
  )
  score0_true <- observed_score(risk_untreated, df_cf0, df_cf0$status, metrics)

  score1 <- ip_score_long(
    probabilities = risk_treated,
    data_outcome = df_val_outcome,
    data_long = df_val_long,
    time_horizon = 5,
    treatment_formula = A ~ (A_lag_1 * L),
    treatment_of_interest = rep(1, 5),
    null_model = FALSE,
    metrics = metrics
  )
  score1_true <- observed_score(risk_treated, df_cf1, df_cf1$status, metrics)

  # expect_equal(score0$score, score0_true$score)
  # expect_equal(score1$score, score1_true$score)
  expect_equal(score0$score, score0_true$score, tolerance = 0.01)
  expect_equal(score1$score, score1_true$score, tolerance = 0.01)

})
