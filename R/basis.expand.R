#
# Perform a basis expansion on a n-vector 'x' of size nbasis.
#
basis.expand <- function(x, nbasis,
                         basis.type = c("poly", "trig", "wave")) {

  if (basis.type[1] == "poly") {
    # POLYNOMIAL basis expansion

    nbasis <- check.nbasis.poly(x, nbasis)

    # generate basis exansion PSI of order nbasis
    PSI <- outer(x, 1:nbasis, "^")

  } else if (basis.type[1] == "trig") {
    # TRIGONOMETRIC basis expansion

    # generate and center basis exansion PSI (of order nbasis)
    # Note: The format of the trig basis expansion is changed from
    # the examples used in the original HierBasis paper.
    # 1 - Remove intercept term since we center the data.
    # 2 - Add odd degree trig functions in addition to the even
    #     degree functions.
    PSI <- outer(x, 1:nbasis, FUN = function(z, nb) {
      ifelse(nb %% 2 == 1,
             cos(nb * pi * z),
             sin((nb - 1) * pi * z))
    })
    warning("Warning in hierbasis2::hierbasis(). Parameter 'basis.type = \"trig\"'
            is unfinished and experimental.")

  } else if (basis.type[1] == "wave") {
    # WAVELET basis expansion

    PSI <- outer(x, 0:(nbasis - 1), FUN = function(z, nb) {
      wavelet(z, nb)
    })
    warning("Warning in hierbasis2::hierbasis(). Parameter 'basis.type = \"wave\"'
            is unfinished and experimental.")

  } else {
    stop("Error in hierbasis2::hierbasis(). Parameter 'basis.type'
         only available for 'poly', 'trig', or 'wave' expansions.")
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
