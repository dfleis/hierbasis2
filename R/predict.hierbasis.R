#' Model Predictions for the Univariate hierbasis model
#'
#' The generic S3 method for predictions for objects of class \code{hierbasis}.
#'
#' @param object A fitted object of class '\code{hierbasis}'.
#' @param new.x An optional vector of x values we wish to fit the fitted
#'              functions at. This should be within the range of
#'              the training data.
#' @param interpolate A logical indicator of if we wish to use
#'                    linear interpolation for estimation of fitted values.
#'                    This becomes useful for high dof when the
#'                    estimation of betas on the original scale becomes unstable.
#' @param ... Not used. Other arguments for predict function.
#'
#' @details
#' This function returns a matrix of  predicted values at the specified
#' values of x given by \code{new.x}. Each column corresponds to a lambda value
#' used for fitting the original model.
#'
#' If \code{new.x == NULL} then this function simply returns
#' the fitted values of the estimated function at the original x values used for
#' model fitting. The predicted values are presented for each lambda values.
#'
#' The function also has an option of making predictions
#' via linear interpolation. If \code{TRUE}, a predicted value is equal to the
#' fitted values if \code{new.x} is an x value used for model fitting. For a
#' value between two x values used for model fitting, this simply returns the
#' value of a linear interpolation of the two fitted values.
#'
#' @return
#' \item{fitted.values}{A matrix with \code{length(new.x)} rows and
#'                      \code{nlam} columns}
#'
#' @export
#'
#' @author Annik Gougeon,
#' David Fleischer (\email{david.fleischer@@mail.mcgill.ca}).
#' @references
#' Haris, A., Shojaie, A. and Simon, N. (2016). Nonparametric Regression with
#' Adaptive Smoothness via a Convex Hierarchical Penalty. Available on request
#' by authors.
#'
#' @seealso The original \code{HierBasis} function, as implemented by
#' Haris et al. (2016) can be found via
#' \url{https://github.com/asadharis/HierBasis/}.
#'
predict.hierbasis <- function(object,
                              new.x       = NULL,
                              lam.idx     = NULL,
                              interpolate = FALSE, ...) {
  if (is.null(lam.idx)) {
    lam.idx <- 1:length(object$lambdas)
  }

  if (is.null(new.x)) {
    # return fitted.response if no new predictors are supplied
    object$fitted.values[,lam.idx]

  } else {
    if (!interpolate) {
      basis.expand.out <- basis.expand(new.x, object$nbasis, object$basis.type)
      new.x.expand     <- basis.expand.out$PSI

      # X %*% beta with the intercept
      fitted <- cbind(1, new.x.expand) %*% coef(object, lam.idx = lam.idx)

    } else {
      # return predicted values
      sapply(lam.idx, FUN = function(i) {
        # obtain curve for a particular value
        yhat.temp <- object$fitted.values[, i]
        # return predictions
        approx(x = object$x, y = yhat.temp, xout = new.x)$y
      })
    }

  }
}
