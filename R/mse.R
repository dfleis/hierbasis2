mse <- function(z, zhat) {
  sum(z - zhat)^2/length(z)
}

