rmse <- function(z, zhat) {
  sqrt(sum(z - zhat)^2/length(z))
}

