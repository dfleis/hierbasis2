#
# Initialize cross-validation data structures with a
# first run of a hierbasis/additivehierbasis model.
#
cv.init <- function(x, y, lambdas, ...) {
  n <- length(y)

  if (!is.null(lambdas)) {
    nlam          <- length(lambdas)
    max.lambda    <- max(lambdas)
    lam.min.ratio <- min(lambdas)/max.lambda

    if (class(x) == "matrix") {
      mod.init <- additivehierbasis(x, y,
                                    max.lambda    = max.lambda,
                                    lam.min.ratio = lam.min.ratio,
                                    nlam          = nlam, ...)
    } else {
      mod.init <- hierbasis(x, y,
                            max.lambda    = max.lambda,
                            lam.min.ratio = lam.min.ratio,
                            nlam          = nlam, ...)
    }
  } else {
    if (class(x) == "matrix") {
      mod.init <- additivehierbasis(x, y, ...)
    } else {
      mod.init <- hierbasis(x, y, ...)
    }

    lambdas       <- mod.init$lambdas
    nlam          <- length(lambdas)
    max.lambda    <- max(lambdas)
    lam.min.ratio <- min(lambdas)/max.lambda
  }

  out <- list()
  out$model.fit     <- mod.init
  out$lambdas       <- lambdas
  out$nlam          <- nlam
  out$max.lambda    <- max.lambda
  out$lam.min.ratio <- lam.min.ratio
  out
}
