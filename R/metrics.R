cf_metric <- function(metric, ...) {

  # Handle calplotX pattern
  if (grepl("^calplot[0-9]+$", metric)) {
    x <- as.numeric(sub("^calplot", "", metric))
    return(cf_calplot(..., n = x))
  }

  switch (metric,
    "brier" = cf_brier(...),
    "scaled_brier" = cf_brier_scaled(...),
    "brier_scaled" = cf_brier_scaled(...),
    "auc" = cf_auc(...),
    "oeratio" = cf_oeratio(...),
    "calplot" = cf_calplot(...),
    stop("Performance metric ", metric, " not implemented")
  )
}

# Brier

cf_brier <- function(obs_outcome, cf_pred, pseudo_i, ipw, ...) {

  # cf brier score is brier score of weighted pseudopop
  1 / sum(ipw[pseudo_i]) *
    sum(
      ( cf_pred[pseudo_i] - obs_outcome[pseudo_i] )^2 * ipw[pseudo_i]
    )
}


cf_brier_scaled <- function(obs_outcome, cf_pred, pseudo_i, ipw, ...) {

  brier <- cf_brier(obs_outcome, cf_pred, pseudo_i, ipw)

  nullpred <- stats::weighted.mean(
    x = obs_outcome[pseudo_i],
    w = ipw[pseudo_i]
  )
  nullpreds <- rep(nullpred, length(obs_outcome))


  brier_null <- cf_brier(obs_outcome, nullpreds, pseudo_i, ipw)

  (1 - brier/brier_null)*100

}

# AUC

cf_auc <- function(obs_outcome, cf_pred, pseudo_i, ipw, ...) {
  obs_outcome <- obs_outcome[pseudo_i]
  cf_pred <- cf_pred[pseudo_i]
  ipw <- ipw[pseudo_i]

  stopifnot(
    "no controls (outcome = 0) in pseudopopulation" = 0 %in% obs_outcome,
    "no cases (outcome = 1) in pseudopopulation" = 1 %in% obs_outcome,
    "nonbinary outcome" = length(unique(obs_outcome)) == 2
  )
  o <- order(cf_pred)
  s <- cf_pred[o]
  y <- obs_outcome[o]      # 0/1
  w <- ipw[o]

  # numeric group id per unique score
  g <- cumsum(c(TRUE, diff(s) != 0))

  # aggregate in C
  w1 <- rowsum(w * y, g, reorder = FALSE)[, 1]
  w0 <- rowsum(w * (1 - y), g, reorder = FALSE)[, 1]

  cum_w0 <- cumsum(w0) - w0

  W1 <- sum(w1)
  W0 <- sum(w0)

  if (W1 == 0 || W0 == 0) return(NA_real_)

  sum(w1 * cum_w0 + 0.5 * w1 * w0) / (W1 * W0)
}

# calibration

# # oe ratio
cf_oeratio <- function(obs_outcome, cf_pred, pseudo_i, ipw, ...) {

  observed <- stats::weighted.mean(
    x = obs_outcome[pseudo_i],
    w = ipw[pseudo_i]
  )

  expected <- mean(cf_pred)

  return(observed/expected)
}

# # calibration plot

cf_calplot <- function(obs_outcome, cf_pred, pseudo_i, ipw, n = 8, ...) {

  # it does not make sense to split the data into more groups than there is
  # unique data. E.g. for the null model, we want one group, not 8.
  n_breaks <- min(n, length(unique(cf_pred)))

  cal <- data.frame(obs_outcome, pseudo_i, cf_pred, ipw)
  cal <- cal[order(cf_pred), ]
  if (n_breaks >= 2) {
    cal$group <- cut(seq_len(nrow(cal)), breaks = n_breaks, labels = F)
  } else {
    cal$group <- 1
  }

  cal$group <- factor(cal$group, levels = seq_len(n_breaks))
  mean_preds <- tapply(
    X = cal$cf_pred,
    INDEX = cal$group,
    FUN = mean
  )

  cal_pseudo <- cal[cal$pseudo_i, ]

  mean_obs <- tapply(
    X = cal_pseudo,
    INDEX = cal_pseudo$group,
    FUN = function(x) stats::weighted.mean(x$obs_outcome, x$ipw)
  )

  calplot <- list(pred = unname(mean_preds), obs = unname(mean_obs))

  calplot
}
