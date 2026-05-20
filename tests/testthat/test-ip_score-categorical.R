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

test_that("ip_score categorical treatments work", {

  n <- 10000

  data <- data.frame(id = 1:n)
  data$L <- rnorm(n, 0.2, 5)
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

  model <- glm(Y ~ trt + P, family = "binomial", data = data)

  ip_score(model, data, Y, trt ~ L, "C")

})
