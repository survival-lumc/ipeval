bootstrap_iteration <- function(ip_object, matchcall, call_env) {
  the_call <- matchcall
  # extract the data and predictions from the call/ip_object
  data <- eval(the_call$data, call_env)
  predictions <- ip_object$predictions

  # sample them
  bs_sample <- sample(nrow(data), size = nrow(data), replace = T)
  bs_data <- data[bs_sample, ]
  bs_predictions <- lapply(predictions, function(x) x[bs_sample])

  # rerun the call, this time with bootstrapped data/predictions
  the_call$data <- quote(bs_data)
  the_call$object <- quote(bs_predictions)
  the_call$bootstrap <- 0
  the_call$null_model <- FALSE
  the_call$strip_ipt_models <- TRUE

  score <- eval(
    the_call,
    envir = list(
      bs_data = bs_data,
      bs_predictions = bs_predictions
    ),
    enclos = call_env
  )$score

  return(score)
}

bootstrap <- function(ip_object, matchcall, call_env, iterations, progress) {
  b <- lapply_progress(
    as.list(1:iterations),
    function(x) {
      if (identical(class(ip_object), "ip_score")) {
        return(bootstrap_iteration(ip_object, matchcall, call_env))
      } else {
        stop("unknown class ", class(ip_object), " found for bootstrapping")
      }
    },
    "bootstrapping",
    progress = progress
  )

  # transpose results
  # (iteration > metric > model) -> (metric > model > iteration)

  # for calibration plot:
  # (iteration > metric > [pred/obs, model]) ->
  # (metric > model > iteration > list(pred = , obs = ))
  transposed <- lapply(ip_object$metrics, function(m) {
    P <- lapply(names(ip_object$predictions), function(p) {
      if (m != "calplot") { # 1 numeric result, simple to combine & transpose
        sapply(b, function(i) i[[m]][[p]])
      } else { # calibration plot, consisting of 2 vectors of preds & obs
        lapply(b, function(i) {
          list(
            pred = i[[m]][["pred", p]],
            obs = i[[m]][["obs", p]]
          )
        })
      }
    })
    names(P) <- names(ip_object$predictions)
    P
  })
  names(transposed) <- ip_object$metrics

  # # summarize
  conf.int <- lapply(ip_object$metrics, function(m) {
    CI <- lapply(names(ip_object$predictions), function(p) {
      if (m != "calplot") {
        return(ci(transposed[[m]][[p]], cover = 0.95))
      } else {
        return(NA)
      }
    })
    names(CI) <- names(ip_object$predictions)
    CI
  })
  names(conf.int) <- ip_object$metrics

  list(
    results = conf.int,
    raw = transposed
  )
}

lapply_progress <- function(x, FUN, task_description, progress = TRUE) {
  # same as lapply, but print a progress indicator
  n <- length(x)

  if (progress == FALSE) {
    result <- lapply(as.list(1:n), FUN)
  } else {
    FUN2 <- function(x, i, n) {
      result <- FUN(x)
      cat("\r", task_description, ": ", i, "/", n, "     ")
      return(result)
    }

    result <- lapply(as.list(1:n), function(i) FUN2(x[[i]], i, n))
    cat("\r")
  }

  return(result)
}

ci <- function(values, cover = 0.95) {
  lower <- (1-cover) / 2
  upper <- 1 - lower
  stats::quantile(values, probs = c(lower, upper))
}
