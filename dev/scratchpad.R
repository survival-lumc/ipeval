library(survival)

n <- 1000000

simulate_data <- function(
    n, seed, A = function() {rbinom(n, 1, plogis(0.5+2*L))}
  ) {
  set.seed(seed)
  L <- rnorm(n)
  trt <- A()
  set.seed(seed + 1)
  P <- rnorm(n)
  Y0 <- rbinom(n, 1, plogis(0.5 + L + 1.25 * P))
  Y1 <- rbinom(n, 1, plogis(0.5 + L + 1.25 * P - 0.9))
  Y <- ifelse(trt == 1, Y1, Y0)
  data <- data.frame(id = 1:n, L, A = trt, P, Y, Y0, Y1)
  return(data)
}

df_dev <- simulate_data(1000, 1)
df_val <- simulate_data(n, 2)
df_cf0 <- simulate_data(n, 2, A = function() rep(0, n))
df_cf1 <- simulate_data(n, 2, A = function() rep(1, n))

trt_model <- glm(A ~ L, family = "binomial", data = df_dev)
propensity_score <- predict(trt_model, type = "response")
df_dev$iptw <- 1 / ifelse(df_dev$A == 1, propensity_score, 1 - propensity_score)

causal_model <- glm(Y ~ A + P, family = "binomial", data = df_dev, weights = iptw)
naive_model <- glm(Y ~ A + P, family = "binomial", data = df_dev)
model_P <- glm(Y ~ P, family = "binomial", data = df_dev)
always0 <- rep(0, n)
always1 <- rep(1,n)
random <- runif(n)

models <- list(
  "naive" = naive_model,
  "causal" = causal_model,
  "modelp" = model_P,
  "always0" = always0,
  "always1" = always1,
  "random" = random
)
metrics <- c("auc", "brier", "oeratio")

score0 <-ip_score(models, df_val, Y, A ~ L, 0, quiet = TRUE, metrics = metrics)
observed0 <- observed_score(models, df_cf0, Y, null_model = TRUE, metrics = metrics)

expect_equal(score0$score, observed0$score, tolerance = 0.01)

score1 <-ip_score(models, df_val, Y, A ~ L, 1, quiet = TRUE, metrics = metrics)
observed1 <- observed_score(models, df_cf1, Y, null_model = TRUE, metrics = metrics)

expect_equal(score1$score, observed1$score, tolerance = 0.01)

