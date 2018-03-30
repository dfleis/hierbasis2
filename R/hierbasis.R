#' Univariate Nonparametric Regression Estimation
#' via a Hierarchical Penalty
#'
#' @details Solves the univariate minimization problem (see Haris et al. 2016)
#' \deqn{argmin_{\beta} (1/2n) || y - \Psi_K \beta ||^2_2 + \lambda\Omega(\beta; m),}
#' where \eqn{\beta} is a real-valued vector of length \eqn{K = } \code{nbasis},
#' \eqn{\Psi_K} is a \eqn{n\times K} real-valued matrix whose columns
#' correspond to the basis expansion of \eqn{x} (ordered in increasing)
#' complexity. The penalty function \eqn{\Omega} is defined as
#' \deqn{\Omega(\beta; m) = \sum^K_{k = 1} w_{k, m} || \Psi_{k:K} \beta_{k:K} ||_2,}
#' for \eqn{\Psi_{k:K}} denotes the submatrix of \eqn{\Psi_K} corresponding to
#' the final columns \eqn{k, k + 1, ..., K}, and \eqn{\beta_{k:K}}
#' the subvector of \eqn{\beta} corresponding to the final elements
#' \eqn{k, k + 1, ..., K}. The weights \eqn{w_{k, m}} in \eqn{\Omega}
#' are given by the m - 1 degree polynomial
#' \deqn{w_{k, m} = k^m - (k - 1)^m},
#' controlling the level of decay on higher order estimates of \eqn{\beta}.
#'
#' @param x A univariate vector representing the predictor.
#' @param y A univariate vector respresenting the response.
#' @param nbasis The number of basis functions to be used in the
#' basis expansion of \code{x}. Default: \code{nbasis = length(y)}.
#' @param max.lambda The largest value of \eqn{\lambda} penalizing
#' the hierarchical penalty. Default: \code{max.lambda = NULL}
#' (computed data-adaptively).
#' @param lam.min.ratio The ratio of the smallest value of \eqn{\lambda}
#' to the largest \eqn{\lambda}. Default: \code{lam.min.ratio = 1e-04}.
#' @param nlam The number \eqn{\lambda}'s to compute the \code{hierbasis}
#' estimators (computed on a base-10 \eqn{log} scale). Default:
#' \code{nlam = 50}.
#' @param m Smoothing parameter controlling the degree of polynomial
#' smoothing/decay performed by the weights in the penalty \eqn{\Omega}.
#' Default: \code{m = 3}.
#' @param type Specifies either Gaussian regression (\code{"gaussian"})
#' with a continuous response vector or binomial regression
#' (\code{"binomial"}) with binary response. Default type is assumed to be
#' Gaussian.
#'
#' @return Returns an object of class \code{hierbasis} with elements
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
hierbasis <- function(x, y,
                      nbasis        = length(y),
                      max.lambda    = NULL,
                      lam.min.ratio = 1e-04,
                      nlam          = 50,
                      m             = 3,
                      type          = c("gaussian", "binomial"))
{
  # To do:
  #   * Implement binomial hierbasis estimator.
  #   * Consider multicategory responses (not in HierBasis).
  #   * Consider renaming 'type' to 'family' in order to be consistent.
  #     with other regression functions implemented in R.
  #   * Add max.iter, tol, max.iter.inner, tol.inner parameters.
  #   * Implementing nonpolynomial basis expansion (wavelet, trig, etc).
  #   * Implementing other weight structures.

  # extract number of observations
  n <- length(y)

  # nbasis.max <- floor(logb(.Machine$double.xmax, base = max(abs(x))))
  # if (nbasis > nbasis.max) {
  #   warning(
  #     paste0("Warning in hierbasis2::hierbasis(): ",
  #             "Basis dimension leads to values too ",
  #             "large for R's floating-point arithmetic. ",
  #             "Setting nbasis = ", nbasis.max, ".")
  #   )
  #   nbasis <- nbasis.max
  # }


  # generate and center basis exansion PSI (of order nbasis)
  PSI <- outer(x, 1:nbasis, "^")
  #PSI.c <- center.fast(PSI)
  #PSIbar <- attr(PSI.c, "scaled:center")
  PSI.c <- scale(PSI, scale = F)
  PSIbar <- attributes(PSI.c)[[2]]

  # center responses
  ybar <- mean(y)
  y.c <- y - ybar

  # compute penalty weights
  w <- (1:nbasis)^m - (0:(nbasis - 1))^m

  # weight matrix (includes tuning parameter lambda)
  # organized columnwise (i.e., rows correspond to a single
  # weight w_k, columns corresponds to a single tuning
  # parameter lambda)
  w.mat <- matrix(rep(w, nlam), ncol = nlam)

  if (is.null(max.lambda)) {
    max.lambda <- NA
  }

  # create empty the hierbasis object to be returned
  out <- list(); class(out) <- "HierBasis"

  if (type[1] == "gaussian") {
    # linear regression analogue

    # compute and extract fits
    hierbasis.fit <- solvehierbasis(PSI.c, y.c, w, w.mat, n,
                                    lam.min.ratio, nlam, max.lambda)
    betahat <- hierbasis.fit$beta

    # compute intercepts for each fitted model
    intercept <- as.vector(ybar - PSIbar %*% betahat)

    # get size of active set per lambda (nonzero coefficients)
    active.set <- colSums((betahat != 0) * 1)

    # compute fitted response
    yhat <- apply(PSI %*% betahat, 1, "+", intercept)

    out$x                     <- x
    out$y                     <- y
    out$m                     <- m
    out$nbasis                <- nbasis
    out$basis.expansion       <- PSI
    out$basis.expansion.means <- PSIbar
    out$ybar                  <- ybar

    out$lambda        <- as.vector(hierbasis.fit$lambda)
    out$intercept     <- intercept
    out$beta          <- betahat
    out$fitted.values <- t(yhat)
    out$n.active      <- active.set
    out$call          <- match.call()

  } else if (type[1] == "binomial") {
    # logistic regression analogue

    stop("Error in hierbasis2::hierbasis(). Parameter 'type = \"binomial\"'
         not yet implemented.")
  } else {
    stop("Error in hierbasis2::hierbasis(). Parameter 'type' only available
         for 'gaussian' or 'binomial'")
  }

  return (out)
}










