data_long <- df_val_long
data_outcome <- df_val_outcome
max_visit_time <- tapply(data_long$visit_time, data_long$id, max)
stopifnot(
  "there are subjects with visit time beyond their survival time" =
    all(max_visit_time <=
          data_outcome$time[match(names(max_visit_time), data_outcome$id)])
)

data_long[1, "visit_time"] <- 10
