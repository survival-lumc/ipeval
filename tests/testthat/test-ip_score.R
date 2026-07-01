# input checks

test_that("wrong input throws sensible errors", {
  n <- 1000
  adminstrative_censor <- 10
  my_data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  my_data$A <- rbinom(n, 1, plogis(0.2 + 0.5*my_data$L))

  my_data$time0 <- simulate_time_to_event(n, 0.04, my_data$L + 0.5*my_data$P)
  my_data$time1 <- simulate_time_to_event(n, 0.04, my_data$L + 0.5*my_data$P - 0.6)
  my_data$time_uncensored <- ifelse(my_data$A == 1, my_data$time1, my_data$time0)
  my_data$status <- ifelse(my_data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  my_data$time <- ifelse(my_data$status == TRUE, my_data$time_uncensored, adminstrative_censor)

  predictions <- runif(n, 0, 1)

  # object ------------------------------------------------------------------
  expect_error(
    ip_score(predictions[1:999], data = my_data, outcome = status, A ~ L, 1),
    "Predictions must be of length nrow"
  )
  expect_error(
    ip_score(lm(status ~ A, my_data), data = my_data, outcome = status, A ~ L, 1),
    "model class lm not supported"
  )

  expect_error(
    ip_score(runif(n, -1, 2), data = my_data, outcome = status, A ~ L, 1),
    "Predictions must be in interval"
  )

  # outcome -----------------------------------------------------------------
  expect_error(
    ip_score(predictions, data = my_data, outcome = B, A ~ L, 1),
    "Outcome B not found"
  )

  expect_error(
    ip_score(predictions, data = my_data, outcome = data$status[1:999], A ~ L, 1),
    "Outcome must be of length"
  )

  expect_error(
    ip_score(predictions, data = my_data, outcome = time, A ~ L, 1),
    "Outcome must be binary"
  )

  expect_error(
    ip_score(predictions, data = my_data, outcome = rep(1, 1000), A ~ L, 1,
            metrics = "auc"),
    "no controls"
  )

  # treatment ---------------------------------------------------------------


  expect_error(
    ip_score(predictions, my_data, status, A + B ~ L, 1),
    "treatment formula must be one variable"
  )

  expect_error(
    ip_score(predictions, my_data, status, ~ L, 1),
    "treatment formula must be one variable"
  )

  expect_error(
    ip_score(predictions, my_data, status, time ~ L, 1),
    "not a factor variable."
  )

  expect_error(
    ip_score(predictions, my_data, status, A ~ L, 2),
    "treatment_of_interest value does not appear in data"
  )

  # other
  expect_error(
    ip_score(predictions, my_data, status, A ~ L, 1,
             iptw = runif(n), bootstrap = 50, bootstrap_progress = FALSE),
    "can't bootstrap"
  )
})



test_that("supplying (list of) model or predictions equivalent", {
  set.seed(1)
  n <- 1000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.5+0.2*data$L))
  data$Y <- rbinom(n, 1, plogis(0.3*data$L + 0.6*data$P - 0.5*data$A))

  model1 <- glm(Y ~ P, family = "binomial", data = data)
  model2 <- glm(Y ~ A + P, family = "binomial", data = data)

  expect_equal(
    ip_score(
      data = data,
      object = model1,
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0,
    ),
    ip_score(
      data = data,
      object = list(model1),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )

  expect_equal(
    ip_score(
      data = data,
      object = model2,
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    ),
    ip_score(
      data = data,
      object = list("model2" = predict_CF(model2, data, "A", 0)),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )

  expect_equal(
    ip_score(
      data = data,
      object = list(aa = model2),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    ),
    ip_score(
      data = data,
      object = list(aa = predict_CF(model2, data, "A", 0)),
      outcome = Y,
      treatment_formula = A ~ L,
      treatment_of_interest = 0
    )
  )
})

test_that("survival outcome at fixed timepoint equivalent to binary outcome", {
  # If we have a binary outcome which is measured for everyone at one fixed
  # timepoint (i.e. t = 10), with no censoring, is specifying outcome = binary
  # equivalent to outcome = Surv(time = 10, binary), time_horizon = 10?

  n <- 10000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)
  data$time <- 10

  model <- suppressWarnings(
    glm(
      Y ~ A + P,
      family = "binomial",
      data = data,
      weights = ipt_weights(data, A ~ L)$weights
    )
  )

  expect_equal(
    ip_score(model, data, Y, A ~ L, 0)$score,
    ip_score(model, data, Surv(time, Y), A ~ L, 0, time_horizon = 10)$score
  )

  expect_equal(
    ip_score(model, data, Surv(time, Y), A ~ L, 0, time_horizon = 10)$ipc$weights,
    rep(1, nrow(data))
  )
})

test_that("iptw/ipcw manual specification equivalent to models", {
  n <- 1000
  horizon <- 10
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  predictions <- runif(n, 0, 1)

  cf_model_iptw <- ip_score(predictions, data, status, A ~ L, 0)
  cf_manual_iptw <- ip_score(predictions, data, status, A ~ 1, 0,
                            iptw = cf_model_iptw$ipt$weights)

  expect_equal(cf_model_iptw$score, cf_manual_iptw$score)
  expect_equal(cf_model_iptw$ipt$weights, cf_manual_iptw$ipt$weights)

  cf_model_ipcw <- ip_score(predictions, data, survival::Surv(time, status), A ~ L, 0,
                           time_horizon = 5)
  cf_manual_ipcw <- ip_score(predictions, data, survival::Surv(time, status), A ~ L, 0,
                            time_horizon = 5, ipcw = cf_model_ipcw$ipc$weights)

  expect_equal(cf_model_ipcw$score, cf_manual_ipcw$score)
  expect_equal(cf_model_ipcw$ipc$weights, cf_manual_ipcw$ipc$weights)

  # bootstrap blocked for iptw vector?
  expect_error(
    cf_manual_iptw <- ip_score(predictions, data, status, A ~ 1, 0,
                               iptw = cf_model_iptw$ipt$weights, bootstrap = 5),
    regexp = "can't bootstrap if iptw are given"
  )


  # test if passing a function for iptw works

  my_iptw_1 <- function(data) {
    rep(1, nrow(data))
  }
  my_iptw_bad <- function(data) {
    rep(1, nrow(data) + 1)
  }
  my_iptw_correct <- function(data) {
    ipt_weights(data, A ~ L, treatment_of_interest = 0)$weights
  }

  expect_equal(
    ip_score(predictions, data, status, A ~ L, 0, iptw = my_iptw_1)$ipt$weights,
    rep(1, nrow(data))
  )


  expect_error(
    ip_score(predictions, data, status, A ~ L, 0, iptw = my_iptw_bad),
    "function specified in iptw did not return a numeric vector of length"
  )

  expect_equal(
    ip_score(predictions, data, status, A ~ L, 0, iptw = my_iptw_correct)$score,
    cf_model_iptw$score
  )
  expect_equal(
    ip_score(predictions, data, status, A ~ L, 0, iptw = my_iptw_correct)$ipt$weights,
    cf_model_iptw$ipt$weights
  )

  set.seed(42)
  bootstrap_manual_fct <-
    ip_score(predictions, data, status, A ~ L, 0, iptw = my_iptw_correct,
             bootstrap = 5)
  set.seed(42)
  bootstrap_normal <- ip_score(predictions, data, status, A ~ L, 0, bootstrap = 5)

  expect_equal(bootstrap_manual_fct$bootstrap, bootstrap_normal$bootstrap)

})

# metrics

test_that("ip_score metrics equal to unobserved CF metrics, binary outcome", {

  set.seed(1)
  n <- 100000
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

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = Y,
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    null_model = FALSE,
    metrics = c("auc", "brier", "oeratio", "scaled_brier")
  )

  Y0_predicted <- predict_CF(model, data, "A", 0)
  score <- riskRegression::Score(
    list(Y0_predicted),
    formula = Y0 ~ 1,
    data = data,
    null.model = T
  )

  brier_model <- score$Brier$score[[2, "Brier"]]
  brier_null <- score$Brier$score[[1, "Brier"]]
  scaled_brier <- (1 - brier_model/brier_null)*100


  score$oe <- mean(data$Y0)/mean(Y0_predicted)


  expect_equal(unname(ip_score$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score$score$brier), score$Brier$score$Brier[[2]], tolerance = 0.01)
  expect_equal(unname(ip_score$score$oeratio), score$oe, tolerance = 0.01)
  expect_equal(unname(ip_score$score$scaled_brier), scaled_brier, tolerance = 0.02)
})

test_that("ip_score metrics equal to unobserved CF metrics, surv, uncensored", {
  set.seed(1)
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))
  data$status <- 1 # no censoring, so status is always 1 at the end

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time <- ifelse(data$A == 1, data$time1, data$time0)

  summary(data$time0)
  summary(data$time1)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data,
    weights = ipt_weights(data, A ~ L)$weights
  )

  horizon <- 10

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(ip_score$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("ip_score metrics equal to unobserved CF metrics, surv, censor at T", {
  set.seed(1)
  horizon <- 9.999
  adminstrative_censor <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)

  summary(data$time0)
  summary(data$time1)

  data$status <- ifelse(data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE, data$time_uncensored, adminstrative_censor)

  data$status_uncensored <- 1

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  ) # naive model, i.e. not adjusting for confounding

  ip_score_km <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  ip_score_cox <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  expect_equal(
    ip_score_km$ipc$weights,
    rep(1, n)
  )
  expect_equal(
    ip_score_cox$ipc$weights,
    rep(1, n)
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(ip_score_km$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score_km$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score_km$score$oeratio), score$oe, tolerance = 0.01)

  expect_equal(unname(ip_score_cox$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score_cox$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score_cox$score$oeratio), score$oe, tolerance = 0.01)
})

test_that("ip_score metrics equal to unobserved CF metrics, surv, KM censor", {
  set.seed(1)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  summary(data$time0)
  summary(data$time1)
  summary(data$censortime)

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(ip_score$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score$score$oeratio), score$oe, tolerance = 0.01)

  # also try treatment == 1 for good measure
  ip_score1 <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 1,
    time_horizon = horizon,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = FALSE
  )

  time1_predicted <- predict_CF(model, data, "A", 1, horizon)
  score1 <- riskRegression::Score(
    list(time1_predicted),
    formula = Hist(time1, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score1$oe <- mean(data$time1 <= horizon)/mean(time1_predicted)

  expect_equal(unname(ip_score1$score$auc), score1$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score1$score$brier), score1$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score1$score$oeratio), score1$oe, tolerance = 0.01)
})

test_that("ip_score metrics equal to unobserved CF metrics, surv, cox censor", {
  set.seed(1)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  summary(data$time0)
  summary(data$time1)
  summary(data$censortime)

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = horizon,
    cens_model = "cox",
    cens_formula = ~ L + P + A,
    null_model = FALSE
  )

  time0_predicted <- predict_CF(model, data, "A", 0, horizon)
  score <- riskRegression::Score(
    list(time0_predicted),
    formula = Hist(time0, status_uncensored) ~ 1,
    data = data,
    null.model = F,
    times = horizon
  )
  score$oe <- mean(data$time0 <= horizon)/mean(time0_predicted)

  expect_equal(unname(ip_score$score$auc), score$AUC$score$AUC, tolerance = 0.01)
  expect_equal(unname(ip_score$score$brier), score$Brier$score$Brier, tolerance = 0.01)
  expect_equal(unname(ip_score$score$oeratio), score$oe, tolerance = 0.01)
})


# minor bootstrap tests
test_that("results are in between lower & upper bootstrap", {
  set.seed(1)
  n <- 1000
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
      data = data
    )
  )

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = Y,
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    bootstrap = 200,
    null_model = FALSE,
    bootstrap_progress = FALSE
  )

  expect_true(
    ip_score$score$auc > ip_score$bootstrap$results$auc[[1]][1] &
      ip_score$score$auc < ip_score$bootstrap$results$auc[[1]][2]
  )
  expect_true(
    ip_score$score$brier > ip_score$bootstrap$results$brier[[1]][1] &
      ip_score$score$brier < ip_score$bootstrap$results$brier[[1]][2]
  )
  expect_true(
    ip_score$score$oeratio > ip_score$bootstrap$results$oeratio[[1]][1] &
      ip_score$score$oeratio < ip_score$bootstrap$results$oeratio[[1]][2]
  )
})

test_that("results are in between lower & upper bootstrap, surv, cox censor", {
  set.seed(1)
  horizon <- 10
  n <- 10000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )
  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    cens_model = "cox",
    time_horizon = horizon,
    cens_formula = ~ L + P + A,
    bootstrap = 100,
    bootstrap_progress = FALSE,
    null_model = FALSE
  )

  expect_true(
    ip_score$score$auc > ip_score$bootstrap$results$auc[[1]][1] &
      ip_score$score$auc < ip_score$bootstrap$results$auc[[1]][2]
  )
  expect_true(
    ip_score$score$brier > ip_score$bootstrap$results$brier[[1]][1] &
      ip_score$score$brier < ip_score$bootstrap$results$brier[[1]][2]
  )
  expect_true(
    ip_score$score$oeratio > ip_score$bootstrap$results$oeratio[[1]][1] &
      ip_score$score$oeratio < ip_score$bootstrap$results$oeratio[[1]][2]
  )
})

# null model

test_that("null model binary outcome", {
  set.seed(1)
  n <- 100000
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
      data = data
    )
  )

  nullmodel <- lm(Y0 ~ 1, data = data)

  ip_score <- ip_score(model, data, Y, A ~ L, 0,
                     metrics = c("brier", "scaled_brier"), null_model = TRUE)
  expect_equal(
    unname(nullmodel$coefficients[1]),
    ip_score$predictions$`null model`[[1]],
    tolerance = 0.01
  )

  # scaled brier score of null model is 0
  expect_equal(
    ip_score$score$scaled_brier[[1]],
    0
  )

  # scaled brier score of model is correct
  expect_equal(
    (1 - ip_score$score$brier[[2]]/ip_score$score$brier[[1]])*100,
    ip_score$score$scaled_brier[[2]]
  )
})

test_that("null model survival outcome", {
  set.seed(2)
  horizon <- 10
  n <- 100000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$censortime <- simulate_time_to_event(n, 0.04, 0.5*data$L + 0.6*data$P +
                                              0.2*data$A)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)
  data$status_uncensored <- 1

  data$status <- ifelse(data$time_uncensored <= data$censortime, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE,
                      data$time_uncensored,
                      data$censortime)


  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    x = TRUE,
    data = data
  )
  nullmodel <- survival::survfit(survival::Surv(time0, status_uncensored) ~ 1, data)

  ip_score <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    metrics = c(),
    cens_model = "cox",
    cens_formula = ~ L + P + A,
    time_horizon = horizon,
    null_model = TRUE
  )

  expect_equal(
    ip_score$predictions$`null model`[[1]],
    1 - summary(nullmodel, times = horizon)$surv,
    tolerance = 0.01
  )

})

# ties

test_that("ip_score handles horizon on censor time correctly", {
  set.seed(1)

  adminstrative_censor <- 10
  n <- 10000
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2 + 0.5*data$L))

  data$time0 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P)
  data$time1 <- simulate_time_to_event(n, 0.04, data$L + 0.5*data$P - 0.6)
  data$time_uncensored <- ifelse(data$A == 1, data$time1, data$time0)

  data$status <- ifelse(data$time_uncensored <= adminstrative_censor, TRUE, FALSE)
  data$time <- ifelse(data$status == TRUE, data$time_uncensored, adminstrative_censor)

  data$status_uncensored <- 1

  model <- survival::coxph(
    formula = survival::Surv(time, status) ~ P + A,
    data = data
  )

  ip_score_km999 <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 9.999,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  ip_score_km10 <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 10,
    cens_model = "KM",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  ip_score_cox999 <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 9.999,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  ip_score_cox10 <- ip_score(
    data = data,
    object = model,
    outcome = survival::Surv(time, status),
    treatment_formula = A ~ L,
    treatment_of_interest = 0,
    time_horizon = 10,
    cens_model = "cox",
    cens_formula = ~ 1,
    null_model = TRUE
  )

  expect_equal(ip_score_km999$score, ip_score_km10$score)
  expect_equal(ip_score_cox999$score, ip_score_cox10$score)
})





