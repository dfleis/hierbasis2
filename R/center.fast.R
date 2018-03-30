.center.fast <- function(M) {
  M.bar <- colMeans(M)
  M.c <- M - rep(1, nrow(M)) %*% t(M.bar)
  attr(M.c, "scaled:center") <- M.bar
  M.c
}

