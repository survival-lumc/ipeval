# Performance in observed dataset

This function computes the performance of the predictions in the given
data, which may contain a mix of treated and untreated subjects. It
exists only to demonstrate the difference between 'normal' performance
and counterfactual performance. It is not user friendly and should not
be relied on. Consider using riskRegression::Score() as an alternative.

## Usage

``` r
observed_score(
  object,
  data,
  outcome,
  metrics = c("auc", "brier", "scaled_brier", "oeratio", "calplot"),
  time_horizon,
  cens_model = "KM",
  cens_formula = ~1,
  null_model = TRUE,
  ipcw
)
```

## Arguments

- object:

  One of the following three options to be validated:

  - a numeric vector, corresponding to risk predictions

  - a glm model

  - a (named) list, with one or more of the previous 2 options, for
    validating and comparing multiple models at once.

- data:

  A data.frame containing the observed outcome.

- outcome:

  The outcome, to be evaluated within data. This could either be the
  name of a numeric/logical column in data, or a Surv object for
  time-to-event data, e.g. Surv(time, status), if time and status are
  columns in data.

- metrics:

  A character vector specifying which performance metrics to be
  computed. Options are c("auc", "brier", "oeratio", "calplot").

- time_horizon:

  For time to event data, the prediction horizon of interest.

- cens_model:

  Model for estimating inverse probability of censored weights (IPCW).
  Methods currently implemented are Kaplan-Meier ("KM") or Cox ("cox"),
  both applied to the censored times. KM is only supported when the
  right hand side of cens_formula is 1.

- cens_formula:

  Formula for which the r.h.s. determines the censoring probabilities.
  I.e. ~ x1 + x2.

- null_model:

  If TRUE fit a risk prediction model which ignores the covariates and
  predicts the same value for all subjects. For time-to-event outcomes,
  the subjects are 'counterfactually' uncensored (using the IPCW, as
  estimated using the cens_formula, or as given by the ipcw argument).

- ipcw:

  A numeric vector, containing the inverse probability of censor
  weights. These are normally computed using the cens_formula, but they
  can be specified directly via this argument.

## Value

Performance metrics in the observed dataset.

## Examples

``` r
n <- 1000

data <- data.frame(L = rnorm(n), P = rnorm(n))
data$A <- rbinom(n, 1, plogis(data$L))
data$Y <- rbinom(n, 1, plogis(0.1 + 0.5*data$L + 0.7*data$P - 2*data$A))

random <- runif(n, 0, 1)
model <- glm(Y ~ A + P, data = data, family = "binomial")

observed_score(
  object = list("ran" = random, "mod" = model),
  data = data,
  outcome = Y,
  metrics = c("auc", "brier", "oeratio")
)
#> 
#>       model   auc brier oeratio
#>  null model 0.500 0.211   1.000
#>         ran 0.467 0.346   0.612
#>         mod 0.760 0.173   1.000
```
