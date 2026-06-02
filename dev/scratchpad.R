set.seed(5)
n <- 1000
data <- data.frame(id = 1:n)

data <- within(data, {
  # 2 visits at t = 0 and t = 2. Time dependent confounding between A and L
  L0 <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(0.7*L0))
  L1 <- rnorm(n, 0.8*L0 - 0.2*A0)
  A1 <- rbinom(n, 1, plogis(0.7*L1 + 0.6*A0))
  # Quickly generate random survival times as an example
  time <- runif(n, 0, 6)
  status <- rbinom(n, 1, 0.5)
})

# If you have data in this 'wide' form, then usually you can use the convenience
# functions to be able to run ip_score_long without doing a lot of extra
# processing
head(data)
# note that due to our simulation there may be measurements after a subject was
# censored or had event. This will be removed later.

# Data preperation for ip_score_long
data_outcome <- data[, c("id", "time", "status")]
data_long <- wide_to_long(
  df = data,
  baseline_variables = "id",
  wide_variables = list("A" = c("A0", "A1"),
                        "L" = c("L0", "L1")),
  visit_times = c(0,2),
  outcome_times = data$time
)
data_long <- add_lag_terms(data_long, "A")

head(data_long)
# note that the measurements that were generated after survival time have been
# removed by wide_to_long().

# To validate, we need predictions. To get started, here is a model that
# randomly predicts a number between 0.2 and 0.8
random_model <- runif(n, min = 0.2, max = 0.8)

# Our treatment at visit i depends on confounder value L at same visit and,
# if visit = 2, on treatment level at previous visit. Thus, we will estimate our
# IPTW model with A ~ L + A_lag_1 in our long data.

# Estimate performance of our random predictions in a pseudopopulation in
# which everybody was assigned treatment strategy {0, 0}
ip_score_long(
  probabilities = random_model,
  data_outcome = data_outcome,
  data_long = data_long,
  treatment_formula = A ~ L + A_lag_1,
  treatment_of_interest = c(0,0),
  visit_times = c(0,2),
  time_horizon = 4
)

# Performance of our random predictions in a pseudopopulation in
# which everybody was assigned treatment level 1 at visit 1, and whatever
# they would normally have at visit 2.
ip_score_long(
  probabilities = random_model,
  data_outcome = data_outcome,
  data_long = data_long,
  treatment_formula = A ~ L + A_lag_1,
  treatment_of_interest = c(1,NA),
  visit_times = c(0,2),
  time_horizon = 4
)
