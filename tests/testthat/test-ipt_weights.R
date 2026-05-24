test_that("IPW function works under binary point trt/binary outcome", {

  build_data <- function(n) {
    df <- data.frame(id = 1:n)
    df$L <- stats::rnorm(n)
    df$A <- stats::rbinom(n, 1, stats::plogis(0.5 + df$L))
    df$P <- stats::rnorm(n)
    df$Y0 <- stats::rbinom(n, 1, stats::plogis(0.5 + df$L + 1.25 * df$P))
    df$Y1 <- stats::rbinom(n, 1, stats::plogis(0.5 + df$L + 1.25 * df$P - 0.6))
    df$Y <- ifelse(df$A == 1, df$Y1, df$Y0)
    return(df)
  }

  set.seed(123)
  df_dev <- build_data(5000)

  in1 <- df_dev$A == 1
  in0 <- df_dev$A == 0


  expect_equal(
    ipt_weights(df_dev, A ~ L, 0)$weights[in0],
    ipw::ipwpoint(A, family = "binomial", link = "logit",
                  denominator =  ~ L, data = df_dev)$ipw.weights[in0]
  )

  expect_equal(
    ipt_weights(df_dev, A ~ L, 1)$weights[in1],
    ipw::ipwpoint(A, family = "binomial", link = "logit",
                  denominator =  ~ L, data = df_dev)$ipw.weights[in1]
  )

})

