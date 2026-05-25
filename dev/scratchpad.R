n <- 10000



data <- data.frame(id = 1:n)
data$L <- rnorm(n, 0.2, 5)
eta2 <- 0.3 + 1.0*data$L
eta3 <- -0.5 - 0.8*data$L

den <- 1 + exp(eta2) + exp(eta3)

p1 <- 1/den
p2 <- exp(eta2)/den
p3 <- exp(eta3)/den


data$A <- sapply(1:n, function(i) {
  sample(c(0,1,2), size=1,
         prob=c(p1[i], p2[i], p3[i]))
})

model <- nnet::multinom(A ~ L, data = data)

predict(model, type = "probs")



class(model)

predict.multinom()


mean(predict(model) == 0)



data$Y0 <- rbinom(n, 1, plogis(0.5 + data$L))
data$Y1 <- rbinom(n, 1, plogis(0.5 + data$L - 0.9))
data$Y2 <- rbinom(n, 1, plogis(0.5 + data$L - 0.5))

data
data$Y <- ifelse(
  test = data$A == 0,
  yes = data$Y0,
  no = ifelse(
    test = data$A == 1,
    yes = data$Y1,
    no = data$Y2
  )
)

