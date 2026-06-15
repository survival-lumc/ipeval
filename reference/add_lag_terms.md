# Add lagged versions of a longitudinal variable

Creates one or more lagged versions of a variable in long-format data.
Lagged values are computed separately within each subject and added as
new columns. The data must be sorted by id and visit time before calling
this function.

## Usage

``` r
add_lag_terms(df, var, lag = 1, fill = 0)
```

## Arguments

- df:

  A long-format data frame containing an \`id\` column identifying
  subjects.

- var:

  Character string giving the name of the variable to lag.

- lag:

  Integer vector specifying the lag(s) to create.

- fill:

  Value used when a lagged observation is unavailable (e.g. at the first
  visit). Defaults to 0.

## Value

The input data frame with additional columns named \`(var)\_lag\_(lag)\`
containing the requested lagged values.

## Examples

``` r
df <- data.frame(
  id = c(1, 1, 1, 2, 2, 2),
  visit_time = c(0, 1, 2, 0, 1, 2),
  A = c(0, 1, 1, 1, 0, 0)
)

# Add a 1-visit lag of A
add_lag_terms(df, "A")
#>   id visit_time A A_lag_1
#> 1  1          0 0       0
#> 2  1          1 1       0
#> 3  1          2 1       1
#> 4  2          0 1       0
#> 5  2          1 0       1
#> 6  2          2 0       0

# Add both 1- and 2-visit lags
add_lag_terms(df, "A", lag = c(1, 2))
#>   id visit_time A A_lag_1 A_lag_2
#> 1  1          0 0       0       0
#> 2  1          1 1       0       0
#> 3  1          2 1       1       0
#> 4  2          0 1       0       0
#> 5  2          1 0       1       0
#> 6  2          2 0       0       1
```
