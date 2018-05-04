expand.poly <- function(x, nbasis) {
  nbasis#nbasis <- check.nbasis.poly(x, nbasis)

  # generate basis exansion PSI of order 1 to nbasis
  PSI <- outer(x, 1:nbasis, "^")
}
