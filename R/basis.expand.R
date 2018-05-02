#
# Perform a basis expansion on a n-vector 'x' of size nbasis.
#
basis.expand <- function(x, nbasis,
                         basis.type = c("poly", "trig", "wave")) {

  if (basis.type[1] == "poly") {
    # POLYNOMIAL basis expansion
    PSI <- expand.poly(x, nbasis)

  } else if (basis.type[1] == "trig") {
    # TRIGONOMETRIC basis expansion
    PSI <- expand.trig(x, nbasis)

  } else if (basis.type[1] == "wave") {
    # WAVELET basis expansion
    PSI <- expand.wave(x, nbasis)

  } else {
    stop("Parameter 'basis.type' only available for 'poly', 'trig', or 'wave' expansions.")
  }

  PSI.c <- scale(PSI, scale = F)
  PSIbar <- attributes(PSI.c)[[2]]

  out <- list()
  out$x           <- x
  out$PSI         <- PSI
  out$PSI.c       <- PSI.c
  out$PSIbar      <- PSIbar
  out$nbasis      <- nbasis
  out$basis.type  <- basis.type
  out
}
