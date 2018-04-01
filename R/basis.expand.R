#
#
#
basis.expand <- function(x, nbasis,
                         basis.type = c("poly", "trig", "wave")) {
  n <- nrow(x)
  p <- ncol(x)

  if (basis.type[1] == "poly") {
    # POLYNOMIAL basis expansion

    # check if the values of the polynomial basis expansion will
    # exceed R's floating point abilities (i.e. max(x)^nbasis is too large)
    nbasis.max <- floor(logb(.Machine$double.xmax, base = max(abs(x))))
    if (nbasis > nbasis.max) {
      warning(
        paste0("Warning in hierbasis2::hierbasis(): ",
               "Basis dimension implies expansion values too ",
               "large for R's floating-point arithmetic. ",
               "Setting nbasis = ", nbasis.max, ".")
      )
      nbasis <- nbasis.max
    }

    # generate and center basis exansion PSI (of order nbasis)
    PSI <- outer(x, 1:nbasis, "^")
    PSI.c <- scale(PSI, scale = F)
    PSIbar <- attributes(PSI.c)[[2]]

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
    PSI.c <- scale(PSI, scale = F)
    PSIbar <- attributes(PSI.c)[[2]]

    warning("Warning in hierbasis2::hierbasis(). Parameter 'basis.type = \"trig\"'
            is unfinished and experimental.")

  } else if (basis.type[1] == "wave") {
    # WAVELET basis expansion

    PSI <- outer(x, 0:(nbasis - 1), FUN = function(z, nb) {
      wavelet(z, nb)
    })
    PSI.c <- scale(PSI, scale = F)
    PSIbar <- attributes(PSI.c)[[2]]

    warning("Warning in hierbasis2::hierbasis(). Parameter 'basis.type = \"wave\"'
            is unfinished and experimental.")

  } else {
    stop("Error in hierbasis2::hierbasis(). Parameter 'basis.type'
         only available for 'poly', 'trig', or 'wave' expansions.")
  }

  out <- list()
  out$x           <- x
  out$PSI         <- PSI
  out$PSI.c       <- PSI.c
  out$PSIbar      <- PSIbar
  out$nbasis      <- nbasis
  out$basis.type  <- basis.type
  out
}
