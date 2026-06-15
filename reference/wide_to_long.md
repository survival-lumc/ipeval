# Convert wide longitudinal data to long format

Reshapes repeated measurements stored in wide format into a long-format
dataset with one row per subject per visit. Rows corresponding to visits
occurring after the subject's outcome time are removed.

## Usage

``` r
wide_to_long(
  df,
  baseline_variables = c("id"),
  wide_variables,
  visit_times,
  outcome_times
)
```

## Arguments

- df:

  A data frame containing one row per subject.

- baseline_variables:

  Character vector of baseline variables to retain and repeat across
  visits. Must include \`"id"\`.

- wide_variables:

  Named list mapping long-format variable names to the corresponding
  wide-format column names for each visit.

- visit_times:

  Numeric vector giving the visit times corresponding to the repeated
  measurements in \`wide_variables\`.

- outcome_times:

  Numeric vector of outcome or follow-up times, used to remove visits
  occurring after a subject's observed follow-up.

## Value

A data frame in long format containing one row per subject per observed
visit, with columns \`id\`, \`visit_time\`, the specified baseline
variables, and the reshaped longitudinal variables.

## Examples

``` r
data <- data.frame(
  id = 1:3,
  A0 = c(0, 1, 0),
  A1 = c(1, NA, 0),
  L0 = c(0.2, -1.1, 0.5),
  L1 = c(0.8, NA, 0.1),
  time = c(3, 1, 4)
)

wide_to_long(
  df = data,
  baseline_variables = "id",
  wide_variables = list(
    A = c("A0", "A1"),
    L = c("L0", "L1")
  ),
  visit_times = c(0, 2),
  outcome_times = data$time
)
#>   id visit_time A    L
#> 1  1          0 0  0.2
#> 2  1          2 1  0.8
#> 3  2          0 1 -1.1
#> 4  3          0 0  0.5
#> 5  3          2 0  0.1
```
