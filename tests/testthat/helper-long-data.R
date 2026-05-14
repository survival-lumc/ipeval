# provides one function to simulate longitudinal data

generate_long_data_cox <- function(
    n, n_visits = 5, seed = 123,
    censoring = FALSE,
    gamma_0 = -1,gamma_L = 0.5,alpha_0 = -2,alpha_A = -0.5,alpha_L = 0.5,
    alpha_U = 0.5,
    U = function() rnorm(n, 0, 0.1),
    Li = function(i) {
      if (i == 1) {
        rnorm(n, U, 1)
      } else {
        rnorm(n, 0.8 * L[, i - 1] - A[, i - 1] + 0.1 * (i-1) + U)
      }
    },
    Ai = function(i) {
      if (i == 1) {
        rbinom(n, 1, plogis(gamma_0 + gamma_L * L[, 1]))
      } else {
        ifelse(
          A[, i - 1] == 1,
          rep(1, n),
          rbinom(n, 1, plogis(gamma_0 + gamma_L * L[, i]))
        )
      }
    },
    coxLP = function(i) {
      alpha_0 + alpha_A * A[, i] + alpha_L * L[, i] + alpha_U * U
    },
    censorLP = function(i) {
      -2 # constant, so non informative. mean censor time is ~ 5
    }
    ) {

  # TODO: add support for setting visit times other than 0, 1, 2, ...
  # visit 1 is always t = 0 though.
  visit_times <- seq(0, n_visits-1)

  L <- matrix(nrow = n, ncol = n_visits)
  A <- matrix(nrow = n, ncol = n_visits)
  time <- rep(NA, n)
  status <- rep(NA, n)

  set.seed(seed)
  U <- U()

  for (i in 1:n_visits) {
    # I have to reset seed every time, because L, A, coxLP may be deterministic
    # or a random process. if its random, then future random draws are different
    # than if it were deterministic. This makes sure that setting trt
    # deterministically to e.g. 0 will make sure outcomes remain consistent
    set.seed(seed + i*100 + 1)

    # simulate A, L, even when we know a patient already has
    # had an event (and therefore misses subsequent visits). These
    # are set to NA later

    L[, i] <- Li(i)

    set.seed(seed + i*100 + 2)
    A[, i] <- Ai(i)

    set.seed(seed + i*100 + 3)
    new_event_time <- simulate_time_to_event(
      n = n,
      constant_baseline_haz = 1,
      LP = coxLP(i)
    )

    if (censoring == TRUE) {
      set.seed(seed + i*100 + 4)
      new_censor_time <- simulate_time_to_event(
        n = n,
        constant_baseline_haz = 1,
        LP = censorLP(i)
      )
      new_time <- pmin(new_event_time, new_censor_time)
      new_status <- new_event_time < new_censor_time

      # abuse memoryless property
      status <- ifelse(
        is.na(status) & new_time < 1,
        new_status,
        status
      )
      time <- ifelse(
        is.na(time) & new_time < 1,
        i - 1 + new_time,
        time
      )
    } else {
      status <- ifelse(
        is.na(status) & new_event_time < 1,
        rep(1, n),
        status
      )
      time <- ifelse(
        is.na(time) & new_event_time < 1,
        i - 1 + new_event_time,
        time
      )
    }
  }

  if (censoring == FALSE) {
    status <- rep(1, n)
    time <- ifelse(is.na(time), n_visits - 1 + new_event_time, time)
  } else {
    status <- ifelse(is.na(status), new_status, status)
    time <- ifelse(is.na(time), n_visits - 1 + new_time, time)
  }

  # wipe A and L values after events (no visit after event)
  for (i in 1:n_visits) {
    A[, i] <- ifelse(time < visit_times[i], rep(NA, n), A[, i])
    L[, i] <- ifelse(time < visit_times[i], rep(NA, n), L[, i])
  }

  colnames(A) <- paste0("A", 0:(n_visits - 1))
  colnames(L) <- paste0("L", 0:(n_visits - 1))

  data.frame(id = 1:n, time, status, A, L, U)
}

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

fit_long_cox_model <- function(data_long, visit_weights) {

  cox_msm <- survival::coxph(
    survival::Surv(time_start, time_end, status) ~ A + A_lag_1 + A_lag_2 +
      A_lag_3 + A_lag_4 + L0,
    data = data_long,
    weights = visit_weights
  )
}

risk_under_0 <- function(cox_msm, t, L0) {
  cumhaz <- survival::basehaz(cox_msm, centered = FALSE)$hazard
  event.times <- survival::basehaz(cox_msm, centered = FALSE)$time
  cumhaz.fun <- stepfun(event.times, c(0, cumhaz))

  # hazard rate term
  hrt <- function(variable, value) {
    exp(coef(cox_msm)[variable]*value)
  }

  # # risk under never treated
  risk0 <- function(t, L0) {
    1 - exp(-cumhaz.fun(t)*hrt("L0", L0))
  }

  return(risk0(t, L0))
}

risk_under_1 <- function(cox_msm, t, L0) {
  cumhaz <- survival::basehaz(cox_msm, centered = FALSE)$hazard
  event.times <- survival::basehaz(cox_msm, centered = FALSE)$time
  cumhaz.fun <- stepfun(event.times, c(0, cumhaz))

  # hazard rate term
  hrt <- function(variable, value) {
    exp(coef(cox_msm)[variable]*value)
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
  return(risk1(t, L0))

}
