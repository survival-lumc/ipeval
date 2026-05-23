test_that("test if assumption output and graphics is correct", {
  set.seed(123)
  n <- 10000
  df <- data.frame(
    conf1 = rnorm(n, mean = 0),
    conf2 = rnorm(n, mean = 0),
    prog = rnorm(n, mean = 0),
    cens1 = runif(n),
    cens2 = runif(n)
  )


  df$trt <- rbinom(n, 1, plogis(0.2 + 0.5*df$conf1 + 0.3*df$conf2))
  df$out <- rbinom(n, 1, plogis(0.3*df$conf1 + 0.6*df$prog - 0.5*df$trt))

  pred <- glm(out ~ conf1 + prog + trt, data = df)

  # correct output for binary point trt, binary outcome
  ips <- ip_score(pred, df, out, trt ~ conf1 + conf2, 1, metrics = c("brier"))
  output <- paste0(capture_output_lines(print(ips)), collapse = "")
  expect_equal(output, "Estimation of the performance of the prediction model in a pseudopopulation where everyone's treatment trt was set to 1.The pseudopopulation ($pseudopop) is constructed from 5469 (54.7%) subjects who indeed received treatment level 1. These subjects are reweighted to represent the full target population under a hypothetical intervention in which everyone received this treatment level.The following assumptions must be satisfied for correct inference:Causal assumptions:- Conditional exchangeability: after adjustment for the covariates used to construct the inverse probability of treatment weights (IPTW), i.e., {conf1, conf2}, there are no unmeasured confounders of treatment assignment and outcome.- Conditional positivity: the probability of receiving treatment level 1 should be greated than zero for each combination of the variables {conf1, conf2} that is observed in the full population. The distribution of IPT-weights can be assessed with $ipt$weights[$pseudopop$ids].- Consistency: the observed outcome under the received treatment equals the potential outcome under that treatment. This includes the assumption of no interference between subjects.Modeling assumptions:- Correctly specified propensity model. Estimated treatment model islogit(trt) = 0.21 + 0.49*conf1 + 0.31*conf2. See also $ipt$model.Performance estimates:      model brier null model 0.239    model.1 0.216")

  # correct output for binary point trt, surv outcome
  ips <- ip_score(pred, df, out, trt ~ conf1 + conf2,time_horizon = todo, 1, metrics = c("brier"))
  output <- paste0(capture_output_lines(print(ips)), collapse = "")
  expect_equal(output, "Estimation of the performance of the prediction model in a pseudopopulation where everyone's treatment trt was set to 1.The pseudopopulation ($pseudopop) is constructed from 5469 (54.7%) subjects who indeed received treatment level 1. These subjects are reweighted to represent the full target population under a hypothetical intervention in which everyone received this treatment level.The following assumptions must be satisfied for correct inference:Causal assumptions:- Conditional exchangeability: after adjustment for the covariates used to construct the inverse probability of treatment weights (IPTW), i.e., {conf1, conf2}, there are no unmeasured confounders of treatment assignment and outcome.- Conditional positivity: the probability of receiving treatment level 1 should be greated than zero for each combination of the variables {conf1, conf2} that is observed in the full population. The distribution of IPT-weights can be assessed with $ipt$weights[$pseudopop$ids].- Consistency: the observed outcome under the received treatment equals the potential outcome under that treatment. This includes the assumption of no interference between subjects.Modeling assumptions:- Correctly specified propensity model. Estimated treatment model islogit(trt) = 0.21 + 0.49*conf1 + 0.31*conf2. See also $ipt$model.Performance estimates:      model brier null model 0.239    model.1 0.216")


  cal <- ip_score(pred, df, out, trt ~ conf, 1, metrics = c("calplot"))

  vdiffr::expect_doppelganger(
    "calibrationplot correct",
    fig = function() plot(cal)
  )
})
