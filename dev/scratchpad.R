df_dev <- generate_long_data_cox(1000, seed = 2)
df_dev_long <- make_dev_long(df_dev)
iptw <- ipt_weights(df_dev_long, A ~ L * A_lag_1)$weights

n <- 500000

coxmsm <- fit_long_cox_model(data_long = df_dev_long, iptw)

df_val <- generate_long_data_cox(n, seed = 3, nonadministrative_censoring = TRUE)
df_cf0 <- generate_long_data_cox(n, seed = 3, nonadministrative_censoring = FALSE,
                                 Af = ~ function() rep(0, n))
df_cf1 <- generate_long_data_cox(n, seed = 3, nonadministrative_censoring = FALSE,
                                 Af = ~ function() rep(1, n))

risk_0 <- risk_under_0(coxmsm, 5, df_val$L0)
risk_1 <- risk_under_1(coxmsm, 5, df_val$L0)

models <- list(risk_0, risk_1)

df_val_outcome <- df_val[, c("id", "time", "status")]
df_val_long <- wide_to_long(
  df_val, "id", list(A = paste0("A", 0:4), L = paste0("L", 0:4)),
  0:4, df_val$time)

df_val_long <- add_lag_terms(df_val_long, "A")

metrics <- c("auc", "brier", "oeratio")


score0_true <- observed_score(models, df_cf0, Surv(time, status),
                              null_model = TRUE, time_horizon = 5,
                              metrics = metrics)

score1_true <- observed_score(models, df_cf1, Surv(time, status),
                              null_model = TRUE, time_horizon = 5,
                              metrics = metrics)




df_val_informativecensor <-
  generate_long_data_cox(n = n, seed = 3, nonadministrative_censoring = TRUE,
                         censorLP = ~ function() -2 + 10*U)


df_val_inf_outcome <- df_val_informativecensor[, c("id", "time", "status")]
df_val_inf_long <- wide_to_long(
  df_val_informativecensor, c("id", "U"), list(A = paste0("A", 0:4), L = paste0("L", 0:4)),
  0:4, df_val_informativecensor$time)

df_val_inf_long <- add_lag_terms(df_val_inf_long, "A")

score0cox_informative <- ip_score_long(
  probabilities = models,
  data_outcome = df_val_inf_outcome,
  data_long = df_val_inf_long,
  visit_times = 0:4,
  time_horizon = 5,
  treatment_formula = A ~ (A_lag_1 * L),
  treatment_of_interest = "never",
  metrics = metrics,
  cens_model = "cox",
  cens_formula = ~ U
)
print_censor_model(score0cox_informative$ipc$model)



# -------------------------------------------------------------------------

survintervals <- survival::survSplit(
  formula = Surv(time, status) ~ .,
  data = df_val_inf_outcome,
  cut = 0:4)

data_combined <- cbind(survintervals, df_val_inf_long)
data_combined$censored <- with(
  data_combined,
  status == 0 & !duplicated(id, fromLast = TRUE))

censorform <- Surv(tstart, time, censored) ~ U
cens_model <- survival::coxph(censorform, data_combined)


bh <- survival::basehaz(cens_model, centered = FALSE)
cumhaz.fun <- stepfun(bh$time, c(0, bh$hazard))

data_combined$cumhaz_start <- cumhaz.fun(data_combined$tstart)
data_combined$cumhaz_end <- cumhaz.fun(pmin(data_combined$time, 5)) # time horizon
data_combined$lp <- predict(cens_model, newdata = data_combined, type = "lp")
data_combined$contribution <- with(data_combined, -(cumhaz_end - cumhaz_start) * exp(lp))

cumhaz <- tapply(data_combined$contribution, data_combined$id, FUN = sum)
prob_uncensor <- exp(cumhaz)

weight <- ifelse(
  df_val_inf_outcome[, "status"] == FALSE & df_val_inf_outcome[, "time"] < 5,
  0,
  1/prob_uncensor)


# -------------------------------------------------------------------------
expect_equal(weight[1:10000], as.vector(score0cox_informative$ipc$weights)[1:10000])

data_combined[, ]
weight[2]
df_val_inf_outcome[2, ]
