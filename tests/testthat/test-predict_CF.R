test_that("predict_CF works for categorical trt", {
  # I wasnt too sure whether the predict CF fct works when handling categorical
  # variables
  set.seed(123)

  n <- 1000

  data <- data.frame(id = 1:n)
  data$L <- rnorm(n, 0.2, 2)
  data$P <- rnorm(n)
  eta2 <- 0.3 + 1.0*data$L
  eta3 <- -0.5 - 0.8*data$L

  den <- 1 + exp(eta2) + exp(eta3)

  p1 <- 1/den
  p2 <- exp(eta2)/den
  p3 <- exp(eta3)/den


  data$trt <- factor(sapply(1:n, function(i) {
    sample(c("A","B","C"), size=1,
           prob=c(p1[i], p2[i], p3[i]))
  }))

  data$YA <- rbinom(n, 1, plogis(2 + data$P - data$L))
  data$YB <- rbinom(n, 1, plogis(2 + data$P - data$L - 0.4))
  data$YC <- rbinom(n, 1, plogis(2 + data$P - data$L - 0.6))
  data$Y <- with(data, ifelse(trt == "A", YA, ifelse(trt == "B", YB, YC)))

  model1 <- glm(Y ~ trt + P, family = "binomial", data = data)
  model2 <- glm(Y ~ P + L, family = "binomial", data = data)

  model1cfa <- predict_CF(model1, data, "trt", "A")
  model1manuala <- plogis(model1$coefficients[[1]] + model1$coefficients[[4]]*data$P)

  model1cfb <- predict_CF(model1, data, "trt", "B")
  model1manualb <- plogis(model1$coefficients[[1]] + model1$coefficients[[2]] + model1$coefficients[[4]]*data$P)

  model1cfc <- predict_CF(model1, data, "trt", "C")
  model1manualc <- plogis(model1$coefficients[[1]] + model1$coefficients[[3]] + model1$coefficients[[4]]*data$P)

  expect_equal(model1cfa, model1manuala)
  expect_equal(model1cfb, model1manualb)
  expect_equal(model1cfc, model1manualc)

  # does the fct work when trt is not part of the model?
  # if it fails, then it likely doesnt work for binary trt either (for which
  # there is no test)
  expect_equal(
    predict_CF(model2, data, "trt", "A"),
    as.vector(predict(model2, type = "response"))
  )

})
