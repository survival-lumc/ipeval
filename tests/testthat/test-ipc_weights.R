test_that("ipc KM simple correct", {
  data <- data.frame(
    my_time = c(1,2,3),
    my_status = c(1,0,1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 0.5)$weights,
    c(1,1,1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 1.5)$weights,
    c(1,1,1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 2)$weights,
    c(1,1,1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 2.5)$weights,
    c(1,NA,2)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 3.5)$weights,
    c(1,NA,2)
  )

  data <- data.frame(
    my_time = c(1,2,3),
    my_status = c(0,1,0)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 0.5)$weights,
    c(1,1,1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 1.5)$weights,
    c(NA,1.5,1.5)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 2.5)$weights,
    c(NA,1.5,1.5)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(my_time, my_status) ~ 1,
                type = "KM", time_horizon = 3.5)$weights,
    c(NA,1.5,NA)
  )

})

# tests with some survival ties (one has outcome, one is censor at same t)
# outcomes (should) happen before censors
test_that("ipc ties", {
  data <- data.frame(
    time = c(1,2,2,3),
    status = c(1, 0, 1, 1)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 1.9)$weights,
    c(1,1,1,1)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 2.0)$weights,
    c(1,1,1,1)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 2.1)$weights,
    c(1,NA,1,2)
  )

  data <- data.frame(
    time = c(1,2,2,3),
    status = c(1, 1, 0, 1)
  )
  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 1.9)$weights,
    c(1,1,1,1)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 2.0)$weights,
    c(1,1,1,1)
  )

  expect_equal(
    ipc_weights(data, survival::Surv(time, status) ~ 1, type = "KM", time_horizon = 2.1)$weights,
    c(1,1,NA,2)
  )

})

test_that(
  "ipc-weighted population represents uncensored population,
  noninformative censoring", {

    horizon <- 20
    build_data <- function(n) {
      df <- data.frame(id = 1:n)
      df$L1 <- stats::rnorm(n)
      df$L2 <- stats::rbinom(n, 1, 0.5)
      df$P1 <- stats::rnorm(n)
      df$P2 <- stats::rbinom(n, 1, 0.5)
      df$A <- rep(1, n)

      u <- runif(n)
      uc <- runif(n)

      lambda0 <- 0.04
      lambdac <- 0.03

      df$failuretime <-  -log(u) / (lambda0 * exp(
        0.5*df$L1 + 0.3*df$L2 + 0.2*df$P1 + 0.4*df$P2))

      # non informative censoring
      df$censortime <- -log(uc)/ (lambdac * exp(0.1))

      df$time <- pmin(df$censortime, df$failuretime)

      df$status <- df$failuretime < df$censortime


      df$time_at_horizon <- pmin(df$time, horizon)
      df$status_at_horizon <- ifelse(df$time < horizon, df$status, F)

      df$status_at_horizon_uncensored <- df$failuretime < horizon
      return(df)
    }
    df_dev <- build_data(1000000)
    df_dev$ipc <- ipc_weights(df_dev, survival::Surv(time, status) ~ 1, "KM", horizon)$weights
    model <- survival::coxph(
      survival::Surv(time_at_horizon, status_at_horizon) ~ L1 + L2 + P1 + P2,
      data = df_dev[!is.na(df_dev$ipc), ], weights = ipc)
    expect_equal(
      unname(model$coefficients),
      c(0.5, 0.3, 0.2, 0.4),
      tolerance = 0.01
    )

    df_dev$pred <- predict_CF(model, df_dev, "A", 1, horizon)

    score <- riskRegression::Score(list(df_dev$pred), Hist(failuretime, A) ~ 1,
                                   data = df_dev, times = horizon, null.model = F)
    score_oeratio <- mean(df_dev$pred)/mean(df_dev$status_at_horizon_uncensored)

    uncensored <- df_dev$time >= horizon | df_dev$status == 1

    ipc_metrics <- with(df_dev, {
      brier <- cf_brier(status_at_horizon, pred, uncensored, ipc)
      auc <- cf_auc(status_at_horizon, pred, uncensored, ipc)
      oeratio <- cf_oeratio(status_at_horizon, pred, uncensored, ipc)
      list(brier, auc, oeratio)
    })


    expect_equal(ipc_metrics[[1]], score$Brier$score$Brier, tolerance = 0.01)
    expect_equal(ipc_metrics[[2]], score$AUC$score$AUC, tolerance = 0.01)
    expect_equal(ipc_metrics[[3]], score_oeratio, tolerance = 0.01)
  })

test_that(
  "ipc-weighted population represents uncensored population,
  cox informative censoring",
  {
    horizon <- 20

    build_data <- function(n) {
      df <- data.frame(id = 1:n)
      df$L1 <- stats::rnorm(n)
      df$L2 <- stats::rbinom(n, 1, 0.5)
      df$P1 <- stats::rnorm(n)
      df$P2 <- stats::rbinom(n, 1, 0.5)
      df$A <- rep(1, n)

      u <- runif(n)
      uc <- runif(n)

      lambda0 <- 0.04
      lambdac <- 0.03

      df$failuretime <-  -log(u) / (lambda0 * exp(
        0.5*df$L1 + 0.3*df$L2 + 0.2*df$P1 + 0.4*df$P2))

      # non informative censoring
      df$censortime <- -log(uc)/ (lambdac * exp(0.4*df$L1 + 0.4*df$P2))

      df$time <- pmin(df$censortime, df$failuretime)
      df$status <- df$failuretime < df$censortime
      df$time_at_horizon <- pmin(df$time, horizon)
      df$status_at_horizon <- ifelse(df$time < horizon, df$status, F)

      df$status_at_horizon_uncensored <- df$failuretime < horizon
      return(df)
    }
    set.seed(1)
    data <- build_data(600000)
    data$ipc <- ipc_weights(data, survival::Surv(time, status) ~ L1 + P2, "cox", horizon)$weights
    model <- survival::coxph(
      survival::Surv(time_at_horizon, status_at_horizon) ~ L1 + L2 + P1 + P2,
      data = data[data$ipc > 0, ], weights = ipc)
    expect_equal(
      unname(model$coefficients),
      c(0.5, 0.3, 0.2, 0.4),
      tolerance = 0.01
    )

    data$pred <- predict_CF(model, data, "A", 1, horizon)

    score <- riskRegression::Score(list(data$pred), Hist(failuretime, A) ~ 1,
                                   data = data, times = horizon, null.model = F)
    score_oeratio <- mean(data$pred)/mean(data$status_at_horizon_uncensored)

    uncensored <- data$time >= horizon | data$status == 1

    ipc_metrics <- with(data, {
      brier <- cf_brier(status_at_horizon, pred, uncensored, ipc)
      auc <- cf_auc(status_at_horizon, pred, uncensored, ipc)
      oeratio <- cf_oeratio(status_at_horizon, pred, uncensored, ipc)
      list(brier, auc, oeratio)
    })

    expect_equal(ipc_metrics[[1]], score$Brier$score$Brier, tolerance = 0.01)
    expect_equal(ipc_metrics[[2]], score$AUC$score$AUC, tolerance = 0.01)
    expect_equal(ipc_metrics[[3]], score_oeratio, tolerance = 0.01)
  }
)
