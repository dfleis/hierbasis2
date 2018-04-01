check.weights <- function(weights, nbasis, m.const) {
  # compute penalty weights
  if (is.null(weights)) {
    w <- (1:nbasis)^m.const - (0:(nbasis - 1))^m.const
  } else {
    # confirm weights has the correct length of nbasis
    if (length(weights) != nbasis) {
      warning("Warning in hierbasis2::hierbasis(): 'weights' vector should be length 'nbasis'.")
    }
    w <- weights
  }

  w
}
