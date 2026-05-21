test_that("ip_score binary treatment but factors", {
  n <- 10000

  set.seed(1)
  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  data$A_factor <- factor(ifelse(data$A == 1, "trt1", "trt0"))

  modelfactor <- glm(Y ~ A_factor + P, data, family = "binomial")
  modelbinary <- glm(Y ~ A + P, data, family = "binomial")

  expect_equal(
    ip_score(modelfactor, data, outcome = Y, A_factor ~ L, "trt0")$score,
    ip_score(modelbinary, data, Y, A ~ L, 0)$score
  )

  expect_equal(
    ip_score(modelfactor, data, Y, A_factor ~ L, "trt1")$score,
    ip_score(modelbinary, data, Y, A ~ L, 1)$score
  )
})

test_that("ip_score categorical treatments binary outcome works", {

  n <- 200000

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

  summary(data[data$trt == "A", "L"])
  summary(data[data$trt == "B", "L"])
  summary(data[data$trt == "C", "L"])

  summary(data$trt)

  data$p1 <- p1
  data$p2 <- p2
  data$p3 <- p3

  data$YA <- rbinom(n, 1, plogis(2 + data$P - data$L))
  data$YB <- rbinom(n, 1, plogis(2 + data$P - data$L - 0.4))
  data$YC <- rbinom(n, 1, plogis(2 + data$P - data$L - 0.6))
  data$Y <- with(data, ifelse(trt == "A", YA, ifelse(trt == "B", YB, YC)))

  mean(data$Y)

  trt_model <- nnet::multinom(trt ~ L, data)

  pred <- predict(trt_model, type = "probs")
  data$estp1 <- unname(pred[, "A"])
  data$estp2 <- unname(pred[, "B"])
  data$estp3 <- unname(pred[, "C"])

  data$weight <- with(
    data,
    ifelse(trt == "A", 1/estp1, ifelse(trt == "B", 1/estp2, 1/estp3))
  )

  mean(data$weight[data$trt == "A"])
  mean(ips$ipt$weights[!is.na(ips$ipt$weights)])

  ips$ipt$weights[!is.na(ips$ipt$weights)][[1]]


  trt_model <- ipt_weights(data, trt ~ L, type = "categorical")
  trt_model$model

  naivemodel <- glm(Y ~ trt + P, family = "binomial", data = data)
  causalmodel <- glm(Y ~ trt + P, family = "binomial", data,
                     weights = trt_model$weights)
  fullmodel <- glm(Y ~ trt + P + L, family = "binomial", data = data)
  predA <- predict_CF(naivemodel, data, "trt", "A")
  predB <- predict_CF(naivemodel, data, "trt", "B")
  predC <- predict_CF(naivemodel, data, "trt", "C")

  naivemodel


  correctA <- data$YA
  correctB <- data$YB
  correctC <- data$YC
  correctnaive <- data$Y
  always0 <- rep(0, n)
  always1 <- rep(1, n)
  random <- runif(n)

  models <- list(predA, predB, predC,
                 correctA, correctB, correctC, correctnaive, always0, always1,
                 random)

  ips <- ip_score(models, data, Y, trt ~ L, "A", null_model = FALSE)
  ips
  observed_score(models, data, YB)



  observed_score(list(naivemodel, causalmodel, fullmodel), data[data$trt != "A",], YB)

  # ip_score(models, data, )

  ips$ipt$weights
  ips$ipt$model





  # observed_score()

  ip_score(list(naivemodel, causalmodel, fullmodel), data[data$trt != "A",],
           outcome = Y, treatment_formula = trt ~ L, "B", null_model = F)
  observed_score(list(naivemodel, causalmodel, fullmodel), data[data$trt != "A",], YB)

})
