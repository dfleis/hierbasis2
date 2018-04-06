check.nbasis.poly <- function(x, nbasis) {
  # check if the values of the polynomial basis expansion will
  # exceed R's floating point abilities (i.e. max(x)^nbasis is too large)
  nbasis.max <- floor(logb(.Machine$double.xmax, base = max(abs(x))))
  if (nbasis > nbasis.max) {
    warning(
      paste0("Warning in check.nbasis.poly(): ",
             "Basis dimension implies expansion values too ",
             "large for R's floating-point arithmetic. ",
             "Setting nbasis = ", nbasis.max, ".")
    )
    nbasis <- nbasis.max
  }
  nbasis
}
