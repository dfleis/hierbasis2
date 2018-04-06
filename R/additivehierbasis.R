#' Multivariate Nonparametric Regression Estimation via
#' a Sparse Hierarchical Penalty
#'
#' The main function for fitting sparse additive models via the additive
#' hierbasis estimator
#'
#' @details Solve the multivariate minimization problem
#' (see Haris at al. (2016) for details) for \eqn{\beta}:
#' \deqn{argmin_{\beta_1,\ldots, \beta_p}
#' (1/(2n)) || y - \sum_j\Psi^{(j)}_K \beta_j ||^2_2 +
#' \lambda^2 (1 - \alpha) \sum_j || \Psi^{(j)}_K \beta_j ||_2 +
#' \lambda \alpha \sum_j \Omega_j( \beta_j; m )  ,}
#' where \eqn{\beta_j} is a vector of length \eqn{K = } \code{nbasis}
#' and the summation is over the index \eqn{j = 1, ..., p} for \eqn{p}
#' covariatres.
#' The penalty function \eqn{\Omega_j(\beta_j; m)} is given by
#' \deqn{\sum_k w_{k, m}\beta_{j, [k:K]},}
#' where \eqn{\beta_{j, [k:K]}} is the vector of coefficients for the
#' \eqn{j}-th predictor and represents \code{beta[k:K]} for the corresponding
#' vector \code{beta}, summing over \eqn{k = 1, ..., K}.
#' Finally, the weights \eqn{w_{k, m}} are given (by default)
#' \deqn{w_{k, m} = k^m - (k - 1)^m,} where \eqn{m} denotes the 'smoothness level'.
#' For details see Haris et al. (2016).
#'
#' @param X An \eqn{n \times p}{n x p} matrix of covariates.
#' @param y A univariate vector respresenting the response.
#' @param nbasis The number of basis functions to be used in the
#'               basis expansion of \code{x}. Default:
#'               \code{nbasis = 10}.
#' @param max.lambda The largest value of \eqn{\lambda} penalizing
#'                   the hierarchical penalty. Default:
#'                   \code{max.lambda = NULL} (computed data-adaptively).
#' @param lam.min.ratio The ratio of the smallest value of \eqn{\lambda}
#'                      to the largest \eqn{\lambda}.
#'                      Default: \code{lam.min.ratio = 1e-04}.
#' @param nlam The number \eqn{\lambda}'s to compute the \code{hierbasis}
#'             estimators (computed on a base-10 \eqn{log} scale). Default:
#'             \code{nlam = 50}.
#' @param beta.mat An initial estimate of the parameter beta, a
#'                 \code{ncol(x)}-by-\code{nbasis} matrix. If NULL
#'                 (default) the inital estimate is set to the zero
#'                 matrix.
#' @param alpha A scalar between 0 and 1 controlling the balance between
#'              the sparsity penalty and the hierarchical penalty.
#'              Default is 0.5.
#' @param m.const Smoothing parameter controlling the degree of polynomial
#'                smoothing/decay performed by the weights in the penalty
#'                \eqn{\Omega}.Default: \code{m.const = 3}. Only relevant if
#'                \code{weights} is left unspecified.
#' @param max.iter Maximum number of iterations for block coordinate descent.
#' @param tol Tolerance/stopping precision for the block coordinate descent algorithm.
#' @param type Specifies either Gaussian regression (\code{"gaussian"})
#'             with a continuous response vector or binomial regression
#'             (\code{"binomial"}) with binary response. Default type is
#'             assumed to be Gaussian.
#' @param weights (New parameter for \code{hierbasis2}) Permits user-specified
#'                weights \eqn{w_{k, m}}. If left unspecified
#'                (\code{weights = NULL}) then the default is set to
#'                \eqn{w_{k, m} = k^m - (k - 1)^m}.
#' @param basis.type (New parameter for \code{hierbasis2}) Desired basis
#'                   functions for the basis expansion of \code{x}. Specifies
#'                   either a polynomial expansion \code{basis.type = "poly"},
#'                   trigonometric expansion \code{basis.type = "trig"}, or
#'                   wavelet expansion \code{basis.type = "wave"}. Default set
#'                   to a polynomial expansion.
#'
#' @return
#' An object of class additivehierbasis with the following elements
#'
#' \item{X, y}{The original \code{X} and \code{y} used as inputs.}
#' \item{beta}{The \code{(nbasis * p)}-by-\code{nlam}-matrix of the estimated
#' linear coefficients \eqn{beta}}.
#' \item{beta.array}{The \code{nbasis * p * nlam}-array of the estimated
#' linear coefficients.}
#' \item{intercept}{The \code{nlam}-vector of estimated intercepts \eqn{\beta_0}.}
#' \item{fitted.values}{The \code{n * nlam} matrix of fitted values.}
#' \item{basis.expansion}{The \code{n * nbasis * p}-array containing the basis
#' expansion of \code{X}. Second-dimension indices correspond to
#' increasing basis function complexity and third-dimension indices correspond
#' to each of the \eqn{p} covariates of \code{X}.}
#' \item{basis.expansion.means}{The \code{nbasis * p}-vector containing the
#' second-dimension-wise means of the basis expansion.}
#' \item{ybar}{Mean of the response vector.}
#' \item{lambdas}{Sequence of tuning parameters \eqn{lambda} used for
#' penalizing the fits.}
#' \item{fitted.values}{The \code{nbasis * nlam}-matrix of fitted responses.}
#' \item{alpha}{The scalar controlling the balance between the sparsity-inducing
#' penalty and the hierarchical-sparsity-inducing penalty.}
#' \item{m.const}{The \code{m.const} value used for defining 'order' of smoothness.}
#' \item{nbasis}{The maximum number of basis functions used for computing
#' the basis expansion of \code{x}.}
#' \item{max.iter}{Maximum number of iterations used for the block coordinate
#' descent algorithm.}
#' \item{tol}{Tolerance/stopping precision used for block coordinate descent.}
#' \item{weights}{The weights used for smoothing/penalizing the basis functions.}
#' \item{active}{The size of the active set (number of nonzero \eqn{\beta})
#' per tuning parameter \eqn{\lambda}.}
#' \item{active.mat}{The size of the active set per predictor (rowwise),
#' per tuning parameter \eqn{\lambda} (columnwise).}
#' \item{type}{The specified family \code{'gaussian'} or \code{'binomial'}.}
#' \item{basis.type}{Specified basis expansion family, polynomial,
#' trigonometric, or wavelet.}
#'
#' @author Annik Gougeon,
#' David Fleischer (\email{david.fleischer@@mail.mcgill.ca}).
#' @references
#' Haris, A., Shojaie, A. and Simon, N. (2016). Nonparametric Regression with
#' Adaptive Smoothness via a Convex Hierarchical Penalty. Available on request
#' by authors.
#'
#' @seealso The original \code{AdditiveHierBasis} function, as implemented by
#' Haris et al. (2016) can be found via
#' \url{https://github.com/asadharis/HierBasis/}.
#'
additivehierbasis <- function(X, y,
                              nbasis        = 10,
                              max.lambda    = NULL,
                              lam.min.ratio = 1e-04,
                              nlam          = 50,
                              beta.mat      = NULL,
                              alpha         = NULL,
                              m.const       = 3,
                              max.iter      = 1e2,
                              tol           = 1e-04,
                              type          = c("gaussian", "binomial"),
                              weights       = NULL,
                              basis.type    = c("poly", "trig", "wave")) {

  # extract sample size, number of covariates, and basis size
  n <- length(y)
  p <- ncol(X)
  K <- nbasis

  if (is.null(beta.mat)) {
    # initialize a matrix of estimates for beta_j
    beta.mat <- matrix(0, ncol = p, nrow = K)
  }

  # compute penalty weights
  w <- check.weights(weights, nbasis, m.const)
  # place into matrix where columns will correspond to varying lambdas
  w.mat <- matrix(rep(w, nlam), ncol = nlam)

  basis.expand.out <- basis.expand.multivar(X, K, basis.type)
  # Note: In the original AdditiveHierBasis code the uncentered
  # design matrix is passed if the "type = 'binomial'", while
  # the centered matrix is used if "type = 'gaussian'"
  nbasis      <- basis.expand.out$nbasis
  PSI.array   <- basis.expand.out$PSI.array
  PSI.c.array <- basis.expand.out$PSI.c
  PSIbar      <- basis.expand.out$PSIbar

  ybar <- mean(y)

  if (is.null(max.lambda)) max.lambda  <- NA
  if (is.null(alpha))      alpha       <- NA

  beta_is_zero <- all(beta.mat == 0)

  if (type[1] == "gaussian") {
    mod <- FitAdditive(y - mean(y),
                       w.mat, w,
                       PSI.c.array, beta.mat,
                       max_lambda    = max.lambda,
                       lam_min_ratio = lam.min.ratio,
                       alpha         = alpha,
                       tol           = tol,
                       p, K, n, nlam,
                       max_iter      = max.iter,
                       beta_is_zero  = beta_is_zero,
                       active_set    = colSums((beta.mat != 0) * 1))
    # betahat is a (nbasis * p)-rows by nlam-columns matrix of fitted
    # coefficients
    betahat       <- mod$beta
    betahat.array <- array(betahat, dim = c(nbasis, p, nlam))

    # compute intercepts
    intercept <- as.vector(ybar - (as.vector(PSIbar) %*% betahat))

    # fit response vector for each lambda
    yhats <- Matrix::crossprod(apply(PSI.c.array, 1, cbind), betahat) + ybar

    # NOTE: What's the best way of computing the active set?
    # Should we just count the total number of nonzero coefficients
    # per tuning parameter lambda? Or should we compute number of nonzero
    # coefficients per predictor per tuning paramter
    active.set <- colSums(as.matrix(betahat) != 0)
    active.set.mat <- apply(betahat.array, 3, function(b) colSums(b != 0))

    # Finally, we return an additivehierbasis object.
    out <- list()
    out$X                     <- X
    out$y                     <- y
    out$beta                  <- betahat
    out$beta.array            <- betahat.array
    out$intercept             <- intercept
    out$fitted.values         <- yhats
    out$basis.expansion       <- PSI.array
    out$basis.expansion.means <- PSIbar
    out$ybar                  <- ybar

    out$lambdas    <- mod$lambdas
    out$alpha      <- alpha
    out$m.const    <- m.const
    out$nbasis     <- nbasis
    out$max.iter   <- max.iter
    out$tol        <- tol
    out$weights    <- w
    out$active     <- active.set
    out$active.mat <- active.set.mat
    out$type       <- type[1]
    out$basis.type <- basis.type[1]
    out$call       <- match.call()

  } else if (type[1] == "binomial") {
    # logistic regression analogue
    stop("Error in hierbasis2::additivehierbasis(). Parameter 'type = \"binomial\"'
         not yet implemented.")

  } else {
    stop("Error in hierbasis2::additivehierbasis(). Parameter 'type' only available
         for 'gaussian' or 'binomial'")
  }

  class(out) <- "additivehierbasis"
  return (out)
}
