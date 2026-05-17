test <- function(n, Af = function(n) {rnorm(n)}, Uf = ~ A + 1) {
  df <- data.frame(id = 1:n)

  A <- Af(n)
  U <- eval(Uf[[2]], envir = list(A = A, n = n))

  df$A <- A
  df$U <- U
  df
}

test(n = 10)
test(n = 10, Af = function(n) {rnorm(n, mean = 10)}, Uf = ~ 1)



df <- generate_long_data_cox(n = 10)
# Li = ~ if (i == 1) { U } else {0.8 * L[, i - 1] - A[, i - 1] + 0.1 * (i-1) + U },
