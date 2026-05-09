n <- 5000
simulate <- function(seed) {
  set.seed(seed)

  data <- data.frame(
    L = rnorm(n, mean = 0),
    P = rnorm(n, mean = 0)
  )
  data$A <- rbinom(n, 1, plogis(0.2+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  data$iptw <- ipt_weights(data, A ~ L)$weights

  y <- data[data$A == 1, "Y"]
  w <- data[data$A == 1, "iptw"]

  estimator <- weighted.mean(y, w)

  truth <- mean(data$Y1)

  return(c("estimator" = estimator, "truth" = truth))
}

results <- lapply_progress(as.list(1:1000), FUN = simulate, task_description = "")
results <- as.data.frame(do.call(rbind, results))
results$bias <- results$estimator - results$truth
t.test(results$bias)
# null model is unbiased



# Q: Is the CF Brier score of the null model, as trained on the CF data,
# an accurate estimator of the true Brier score of the null model, as trained on
# the true data?

simulate2 <- function(seed) {
  set.seed(seed)

  data <- data.frame(
    L = rnorm(n, mean = 1),
    P = rnorm(n, mean = -0.5)
  )
  data$A <- rbinom(n, 1, plogis(-0.3+0.5*data$L))
  data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P))
  data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.3*data$L + 0.4*data$P - 0.4))
  data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

  data$iptw <- ipt_weights(data, A ~ L)$weights

  y <- data[data$A == 1, "Y"]
  w <- data[data$A == 1, "iptw"]

  CFnull <- weighted.mean(y, w)
  truenull <- mean(data$Y1)

  ips <- ip_score(
    list("cfnull" = rep(CFnull, n),
         "truenull" = rep(truenull, n)),
    data, Y, A ~ L, 1, "brier", null_model = FALSE)

  obs <- observed_score(
    list("cfnull" = rep(CFnull, n),
         "truenull" = rep(truenull, n)),
    data, Y1, "brier"
  )
  ips$score$brier
  obs$score$brier

  return(c("cfnullcfdata" = ips$score$brier[[1]],
           "truenullcfdata" = ips$score$brier[[2]],
           "cfnulltruedata" = obs$score$brier[[1]],
           "truenulltruedata" = obs$score$brier[[2]]))
}
# to compare estimate of scaled brier to 'truth', we implicitly use cfnullcfdata
# to estimate truenulltruedata.

results <- lapply_progress(as.list(1:1000), FUN = simulate2, task_description = "")
results <- as.data.frame(do.call(rbind, results))
results$bias <- results$cfnullcfdata - results$truenulltruedata
results$bias2 <- results$cfnullcfdata - results$cfnulltruedata
results$bias3 <- results$cfnullcfdata - results$truenullcfdata
results$bias4 <- results$truenullcfdata - results$truenulltruedata

# due to 'overfitting' or something?
t.test(results$bias) #bias
t.test(results$bias2) # bias
t.test(results$bias3) #bias
t.test(results$bias4) #not

