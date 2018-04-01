predict.cv.hierbasis <- function(cv.object,
                                 new.x       = NULL,
                                 lam.idx     = c("lambda.1se", "lambda.min"),
                                 interpolate = FALSE, ...) {
  nlam <- length(object$lambdas)

  if (is.null(new.x)) {
    # return fitted.response if no new predictors are supplied
    object$fitted.values

  } else {
    if (!interpolate) {
      basis.expand.out <- basis.expand(new.x, object$nbasis, object$basis.type)
      new.x.expand     <- basis.expand.out$PSI

      # X %*% beta without the intercept
      fitted <- new.x.expand %*% object$beta
      # add the intercept
      t(apply(fitted, 1, "+", object$intercept))

    } else {
      # return predicted values
      sapply(1:nlam, FUN = function(i) {
        # obtain curve for a particular value
        yhat.temp <- object$fitted.values[, i]
        # return predictions
        approx(x = object$x, y = yhat.temp, xout = new.x)$y
      })
    }

  }
}
