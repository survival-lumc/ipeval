test_that("test if assumption output and graphics is correct", {
  set.seed(123)
  n <- 10000
  df <- data.frame(
    conf = rnorm(n, mean = 0),
    prog = rnorm(n, mean = 0)
  )
  df$trt <- rbinom(n, 1, plogis(0.2 + 0.5*df$conf))
  df$out <- rbinom(n, 1, plogis(0.3*df$conf + 0.6*df$prog - 0.5*df$trt))

  pred <- glm(out ~ conf + prog + trt, data = df)

  ips <- ip_score(pred, df, out, trt ~ conf, 1, metrics = c("brier"))

  output <- capture_output_lines(print(ips))

  expect_equal(
    paste0(output[1:10], collapse = ""),
    expected <- "Estimation of the performance of the prediction model in a counterfactual (CF) dataset where everyone's treatment trt was set to 1.The following assumptions must be satisfied for correct inference:- Conditional exchangeability requires that the inverse probability of treatment weights are sufficient to adjust for confounding between treatment and outcome.- Conditional positivity (assess $ipt$weights for outliers)- Consistency (including no interference)- Correctly specified propensity model. Estimated treatment model is"
  )
  cal <- ip_score(pred, df, out, trt ~ conf, 1, metrics = c("calplot"))

  vdiffr::expect_doppelganger(
    "calibrationplot correct",
    fig = function() plot(cal)
  )
})
