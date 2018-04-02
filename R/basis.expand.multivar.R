#
# Perform a basis expansion on a (n * p)-matrix X. Each
# column of X will be expanded using nbasis basis functions
# yielding an array output of dimension (n * nbasis * p).
#
basis.expand.multivar <- function(X, nbasis,
                                  basis.type = c("poly", "trig", "wave")) {
  n <- nrow(X)
  p <- ncol(X)

  # each slice of the array has the orthogonal design
  # corresponding to each covariate
  PSI.array <- PSI.c.array <- array(NA, dim = c(n, nbasis, p))

  # initialize an empty matrix to store PSIbar values
  # so that we'll be able to center our data appropriately
  PSIbar <- matrix(NA, ncol = p, nrow = nbasis)

  # generate and center basis exansion PSI over every predictor
  # j = 1, ..., p of order 1 to nbasis
  for (j in 1:p) {
    if (basis.type[1] == "poly") {
      # POLYNOMIAL basis expansion
      PSI.array[,, j] <- expand.poly(X[,j], nbasis)

    } else if (basis.type[1] == "trig") {
      # TRIGONOMETRIC basis expansion
      PSI.array[,, j] <- expand.trig(X[,j], nbasis)

    } else if (basis.type[1] == "wave") {
      # WAVELET basis expansion
      PSI.array[,, j] <- expand.wave(X[,j], nbasis)

    } else {
      stop("Parameter 'basis.type' only available for 'poly', 'trig', or 'wave' expansions.")
    }

    PSIj.c            <- scale(PSI.array[,, j], scale = F)
    PSI.c.array[,, j] <- PSIj.c
    PSIbar[, j]       <- attributes(PSIj.c)[[2]]
  }

  out <- list()
  out$X           <- X
  out$PSI.array   <- PSI.array
  out$PSI.c.array <- PSI.c.array
  out$PSIbar      <- PSIbar
  out$nbasis      <- nbasis
  out$basis.type  <- basis.type
  out
}
