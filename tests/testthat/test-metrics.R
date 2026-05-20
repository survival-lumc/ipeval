test_that(
  "Binary outcome/point trt/simple confounder exactly correct",
  {

    # set-up

    build_data <- function(n) {
      df <- data.frame(id = 1:n)
      df$L <- stats::rbinom(n, 1, 0.5)
      df$A <- stats::rbinom(n, 1, stats::plogis(df$L))
      df$P <- stats::rnorm(n)
      df$Y0 <- stats::rbinom(n, 1, stats::plogis(0.5 + df$L + 1.25 * df$P))
      df$Y1 <- stats::rbinom(n, 1, stats::plogis(0.5 + df$L + 1.25 * df$P - 0.6))
      df$Y <- ifelse(df$A == 1, df$Y1, df$Y0)
      return(df)
    }

    set.seed(1)

    df_toy <- build_data(100)
    # artificially add some ties with flipped outcome to make sure AUC handles
    # ties correctly
    df_flipped <- df_toy[1:10,]
    df_flipped$Y0 <- 1 - df_flipped$Y0
    df_flipped$Y1 <- 1 - df_flipped$Y1
    df_flipped$Y <- 1 - df_flipped$Y
    df_toy <- rbind(df_toy, df_flipped)

    df_toy$ipw <- ipt_weights(df_toy, A ~ L)$weights

    is_whole <- function(x, tol = .Machine$double.eps^0.5) {
      abs(x - round(x)) < tol
    }

    # find smallest a such that x * a is a whole number
    smallest_a_st_xa_int <- function(x) {
      smallest_a_st_xa_int2 <- function(x) {

        for (a in 1:10000) {
          if (is_whole(x*a)) {
            return(a)
          }
        }
        stop("No smallest a <= 10000 found")
      }
      sapply(x, smallest_a_st_xa_int2)
    }

    # least common multiplier
    lcm <- function(x) {
      gcd <- function(a, b) {
        a <- abs(a)
        b <- abs(b)
        while (b != 0) {
          temp <- b
          b <- a %% b
          a <- temp
        }
        a
      }

      lcm2 <- function(a, b) {
        if (a == 0 || b == 0) return(0)
        abs(a * b) / gcd(a, b)
      }
      Reduce(lcm2, x)
    }

    # build an 'unweighted' pseudopop, by repeating each row a whole number
    # proportional to its ipw
    df_pseudo_exact <- df_toy[df_toy$A == 0, ]
    multiply_rows <- lcm(smallest_a_st_xa_int(unique(df_pseudo_exact$ipw)))
    # this needs to be rounded due to rounding errors in R.
    # i.e. 7.0000000001 becomes 7
    df_pseudo_exact$rep <- round(df_pseudo_exact$ipw * multiply_rows)
    df_pseudo_exact <- df_pseudo_exact[rep(
      seq_len(nrow(df_pseudo_exact)),
      times = round(df_pseudo_exact$rep)
    ),]


    # fit a CF model on original toy data
    model <- suppressWarnings(
      glm(Y ~ A + P, family = "binomial",
          data = df_toy, weights = ipw)
    )
    # predict CF outcomed under trt 0
    df_toy$pred <- predict_CF(model, df_toy, "A", 0)
    df_pseudo_exact$pred <- predict_CF(model, df_pseudo_exact, "A", 0)

    score <- riskRegression::Score(list(df_pseudo_exact$pred), Y ~ 1,
                                   data = df_pseudo_exact,
                                   null.model = F, se.fit = F)
    score_oeratio <- mean(df_pseudo_exact$Y)/mean(df_toy$pred)

    cfscore <- with(df_toy, {
      brier <- cf_brier(Y, pred, A == 0, ipw)
      auc <- cf_auc(Y, pred, A == 0, ipw)
      oeratio <- cf_oeratio(Y, pred, A == 0, ipw)
      list(brier, auc, oeratio)
    })

    expect_equal(cfscore[[1]], score$Brier$score$Brier)
    expect_equal(cfscore[[2]], score$AUC$score$AUC)
    expect_equal(cfscore[[3]], score_oeratio)

  })

test_that(
  "Binary outcome/point trt/multiple confounders approx correct",
  {
    build_data <- function(n, shift) {
      df <- data.frame(id = 1:n)
      df$L1 <- stats::rnorm(n)
      df$L2 <- stats::rbinom(n, 1, 0.5) - 0.5
      df$A <- stats::rbinom(n, 1, stats::plogis(shift + df$L1 + df$L2))
      df$P1 <- stats::rnorm(n)
      df$P2 <- stats::rbinom(n, 1, 0.5) - 0.5
      df$Y0 <- stats::rbinom(n, 1, stats::plogis(
        0.5 + df$L1 + df$L2 + 1.25 * df$P1 + df$P2
      ))
      df$Y1 <- stats::rbinom(n, 1, stats::plogis(
        0.5 + df$L1 + df$L2 + 1.25 * df$P1 + df$P2 - 0.6
      ))
      df$Y <- ifelse(df$A == 1, df$Y1, df$Y0)
      return(df)
    }

    set.seed(1)
    df_dev <- build_data(200, shift = 1)
    # causal model
    df_dev$ipw <- ipt_weights(df_dev, A ~ L1 + L2)$weights

    model <- suppressWarnings(
      glm(Y ~ A + P1 + P2, family = "binomial",
          data = df_dev, weights = ipw)
    )

    shift <- -1
    n <- 1000000
    df_val <- build_data(n, shift = shift)
    # we know the ipw formula, so might as well fill that in
    df_val$ipw <- 1/ifelse(
      df_val$A == 0,
      1 - stats::plogis(shift + df_val$L1 + df_val$L2),
      stats::plogis(shift + df_val$L1 + df_val$L2)
    )

    pred0 <- predict_CF(model, df_val, "A", 0)

    # truth vs cf estimation
    truth <- with(df_val, {
      auc <- cf_auc(Y0, pred0, rep(TRUE, n), rep(1,n))
      brier <- cf_brier(Y0, pred0, rep(TRUE, n), rep(1,n))
      oe <- cf_oeratio(Y0, pred0, rep(TRUE, n), rep(1,n))
      calplot <- cf_calplot(Y0, pred0, rep(TRUE, n), rep(1,n))
      list(auc, brier, oe, calplot)
    })

    cf_est <- with(df_val,  {
      auc <- cf_auc(Y, pred0, A == 0, ipw)
      brier <- cf_brier(Y, pred0, A == 0, ipw)
      oe <- cf_oeratio(Y, pred0, A == 0, ipw)
      calplot <- cf_calplot(Y, pred0, A == 0, ipw)
      list(auc, brier, oe, calplot)
    })
    expect_equal(truth, cf_est, tolerance = 0.005)
  })

test_that(
  "binary outcome/point trt trivial auc correct",
  {
    expect_equal(
      cf_auc(
        obs_outcome = c(0,1,0,1),
        cf_pred = c(0.5,0.5,0.5,0.5),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      0.5
    )

    expect_equal(
      cf_auc(
        obs_outcome = c(0,0,1,1),
        cf_pred = c(1,0,1,0),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      0.5
    )
    expect_equal(
      cf_auc(
        obs_outcome = c(0,0,1,1),
        cf_pred = c(0.1,0.2,0.3,0.4),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      1
    )
    expect_equal(
      cf_auc(
        obs_outcome = c(0,0,1,1),
        cf_pred = c(1,1,0,0),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      0
    )

    # todo: expect errors for missing outcomes in trt groups etc.
  }
)

test_that(
  "binary outcome/point trt trivial oe correct",
  {
    expect_equal(
      cf_oeratio(
        obs_outcome = c(0,1,0,1),
        cf_pred = c(0.5,0.5,0.5,0.5),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      1
    )
    expect_equal(
      cf_oeratio(
        obs_outcome = c(0,0,0,0),
        cf_pred = c(0.5,0.5,0.5,0.5),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      0
    )
    expect_equal(
      cf_oeratio(
        obs_outcome = c(1,0,0,0),
        cf_pred = c(0,0.5,0,0),
        pseudo_i = rep(TRUE, 4),
        ipw = c(1,1,1,1)
      ),
      2
    )
  })

test_that(
  "binary outcome/point trt/binary confounder BS analytically correct",
  {
    # this test basically just checks if the correct weighted mean is done
    n <- 1000
    data <- data.frame(
      L = rbinom(n, 1, 0.3)
    )
    data$A <- rbinom(n, 1, plogis(0.2 + 0.6*data$L))
    data$Y0 <- rbinom(n, 1, plogis(0.1 + 0.4*data$L))
    data$Y1 <- rbinom(n, 1, plogis(0.1 + 0.4*data$L - 0.5))
    data$Y <- ifelse(data$A == 1, data$Y1, data$Y0)

    data$predictions <- runif(n, 0, 1)

    p_hat <- with(data, tapply(A, L, mean))
    p_i <- p_hat[as.character(data$L)]
    data$w <- 1/ifelse(data$A == 1, p_i, 1 - p_i)

    expect_equal(
      cf_brier(data$Y, data$predictions, data$A == 0, data$w),
      with(data[data$A == 0,], 1/sum(w) * sum((predictions - Y)^2 * w))
    )
    expect_equal(
      unname(ip_score(data$predictions, data, Y, A ~ L, 0, metrics = "brier",
                     null_model = F)$score$brier),
      with(data[data$A == 0,], 1/sum(w) * sum((predictions - Y)^2 * w))
    )
  }
)

