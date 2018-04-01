#' Model Predictions for the Multivariate additivehierbasis model
#'
#' The generic S3 method for predictions for objects of
#' class \code{additivehierbasis}.
#'
#' @param object A fitted object of class '\code{additivehierbasis}'.
#' @param new.x An optional matrix of values of \code{x} at which to predict.
#' The number of columns of \code{new.x} should be equal to the number of
#' columns of \code{object$x}.
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
#' @return
#' \item{fitted.values}{A matrix of fitted values with \code{nrow(new.x)}
#'                      rows and \code{nlam} columns}
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
predict.additivehierbasis <- function(object, new.X = NULL, ...) {

  # initialize & extract some variables
  if (is.null(new.X)) {
    object$fitted.values
  } else {
    n.new      <- dim(new.X)[1]
    n          <- nrow(object$X)
    p          <- ncol(object$X)
    nlam       <- length(object$lambdas)
    K          <- object$nbasis
    basis.type <- object$basis.type

    # generate design matrices
    basis.expand.out <- basis.expand.multivar(new.X, K, basis.type)
    PSI.array        <- basis.expand.out$PSI.array

    # compute X %*% beta values
    ans <- Matrix::crossprod(apply(PSI.array, 1, cbind), object$beta)

    # add intercept
    ans.out <- t(apply(ans, 1, "+", object$intercept))

    if (object$type[1] == "gaussian") {
      return (ans.out)
    } else {
      1/(1 + exp(-ans.out))
    }
  }
}
