n <- 3000
n_visits <- 5

gamma_0 <- -1
gamma_L <- 0.5

alpha_0 <- -2
alpha_A <- -0.5
alpha_L <- 0.5
alpha_U <- 0.5

U <- rnorm(n, 0, 0.1)

simulate_longitudinal <- function(n, fix_trt = NULL, seed = 123) {
  set.seed(seed)

  A <- matrix(nrow = n, ncol = n_visits)
  L <- matrix(nrow = n, ncol = n_visits)

  time <- rep(NA, n)
  status <- rep(NA, n)

  simulate_A <- function(i, L, fix_trt) {
    if (is.null(fix_trt)) {
      return(rbinom(n, 1, plogis(gamma_0 + gamma_L * L[, i])))
    } else {
      return(rep(fix_trt, n))
    }
  }

  L[, 1] <- rnorm(n, U, 1)
  A[, 1] <- simulate_A(1, L, fix_trt)

  for (i in 2:n_visits) {
    L[, i] <- rnorm(n, 0.8 * L[, i - 1] - A[, i - 1] + 0.1 * (i-1) + U)
    A[, i] <- ifelse(A[, i - 1] == 1, rep(1, n), simulate_A(i, L, fix_trt))
  }

  for (i in 1:n_visits){
    new.t <- ipeval:::simulate_time_to_event(
      n = n,
      constant_baseline_haz = 1,
      LP = alpha_0 + alpha_A * A[, i] + alpha_L * L[, i] + alpha_U * U
    )
    time <- ifelse(is.na(time) & new.t < 1, i - 1 + new.t,time)
  }
  status <- ifelse(is.na(time), 0, 1)
  time <- ifelse(is.na(time), 5, time)

  colnames(A) <- paste0("A", 0:(n_visits - 1))
  colnames(L) <- paste0("L", 0:(n_visits - 1))

  data.frame(id = 1:n, time, status, A, L, U)
}

# prepare into long format for fitting Cox model

make_dev_long <- function(df_dev) {
  df_dev_long <- wide_to_long(
    df_dev,
    baseline_variables = c("id", "time", "status", "L0", "U"),
    wide_variables = list(A = c("A0", "A1", "A2", "A3", "A4"),
                          L = c("L0", "L1", "L2", "L3", "L4")),
    visit_times = c(0,1,2,3,4),
    outcome_times = df_dev$time
  )

  df_dev_long$time_start <- df_dev_long$visit_time
  df_dev_long$time_end <- pmin(df_dev_long$visit_time + 1, df_dev_long$time)
  df_dev_long$status <- ifelse(
    df_dev_long$time_end == df_dev_long$time,
    df_dev_long$status,
    0
  )
  df_dev_long$time <- NULL

  df_dev_long <- add_lag_terms(df_dev_long, "A", 1:4)
  df_dev_long
}


fit_long_models <- function(df_dev_long) {
  ipt_visits <- ipeval:::ipt_weights(df_dev_long, A ~ L * A_lag_1)

  df_dev_long$iptw_visits <- ipt_visits$weights
  df_dev_long$iptw_cumprod <- unlist(tapply(df_dev_long$iptw_visits,
                                            df_dev_long$id, cumprod))

  cox_msm <- coxph(
    Surv(time_start, time_end, status) ~ A + A_lag_1 + A_lag_2 + A_lag_3 + A_lag_4 + L0,
    data = df_dev_long,
    weights = iptw_cumprod
  )

  cumhaz <- basehaz(cox_msm, centered = FALSE)$hazard
  event.times <- basehaz(cox_msm, centered = FALSE)$time
  cumhaz.fun <- stepfun(event.times, c(0, cumhaz))

  # hazard rate term
  hrt <- function(variable, value) {
    exp(coef(cox_msm)[variable]*value)
  }

  # # risk under never treated
  risk0 <- function(t, L0) {
    1 - exp(-cumhaz.fun(t)*hrt("L0", L0))
  }

  # # risk under always treated
  risk1 <- function(t, L0) {
    1 - exp(-(
      cumhaz.fun(min(t, 1))*
        hrt("L0", L0)*hrt("A", 1) +

        (t >= 1) *
        (cumhaz.fun(min(t, 2)) - cumhaz.fun(1))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag_1", 1) +

        (t >= 2) *
        (cumhaz.fun(min(t, 3)) - cumhaz.fun(2))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag_1", 1)*hrt("A_lag_2", 1) +

        (t >= 3) *
        (cumhaz.fun(min(t, 4)) - cumhaz.fun(3))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag_1", 1)*hrt("A_lag_2", 1)*hrt("A_lag_3", 1) +

        (t >= 4) *
        (cumhaz.fun(min(t, 5)) - cumhaz.fun(4))*
        hrt("L0", L0)*hrt("A", 1)*hrt("A_lag_1", 1)*hrt("A_lag_2", 1)*
        hrt("A_lag_3", 1)*hrt("A_lag_4", 1)
    ))
  }
  return(list(risk0, risk1))
}

metrics <- c("oeratio", "auc", "brier", "scaled_brier")

df_dev <- simulate_longitudinal(n = 3000, seed = 1)
df_dev_long <- make_dev_long(df_dev)

models <- tryCatch(
  fit_long_models(df_dev_long),
  error = function(e) NULL
)
if (is.null(models)) {
  cat("model building failed in seed ", seed, "\n")
  return(NULL)
}
risk0 <- models[[1]]
risk1 <- models[[2]]

simulation_run <- function(seed) {
  U <- rnorm(n, 0, 0.1)

  df_val <- simulate_longitudinal(n = 3000, seed = seed + 1)
  predicted0 <- risk0(5, df_val$L0)
  predicted1 <- risk1(5, df_val$L0)

  df_val_outcomes <- df_val[, c("id", "time", "status")]
  df_val_long <- wide_to_long(
    df_val,
    baseline_variables = c("id", "L0"),
    wide_variables = list("A" = paste0("A", 0:4), "L" = paste0("L", 0:4)),
    visit_times = 0:4,
    outcome_times = df_val$time
  )
  df_val_long <- add_lag_terms(df_val_long, "A")

  score0 <- ip_score_long(
    probabilities = predicted0,
    data_outcome = df_val_outcomes,
    data_long = df_val_long,
    time_horizon = 5,
    treatment_formula = A ~ A_lag_1 * L,
    treatment_of_interest = rep(0,5),
    metrics = metrics,
    null_model = TRUE
  )

  nullwm <- score0$predictions$`null model`[[1]]
  nullkm <- score0$predictions$km0[[1]]

  # use same seed as df_val so L0's and predictions coincide. Afterwards
  # treatment is set fixed to never or always to estimate 'truth'.
  df_true0 <- simulate_longitudinal(n = 3000, fix_trt = 0, seed = seed + 1)
  nulltrue <- mean(df_true0$status)

  return(list("biaswm" = nullwm - nulltrue,
              "biaskm" = nullkm - nulltrue))
}

results <- lapply_progress(1:20, simulation_run, "")
results <- do.call(rbind, results)
colMeans(matrix(unlist(results), ncol = 2, byrow = FALSE))

