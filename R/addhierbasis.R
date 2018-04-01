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
#' @param type Specifies either Gaussian regression (\code{"gaussian"})
#'             with a continuous response vector or binomial regression
#'             (\code{"binomial"}) with binary response. Default type is
#'             assumed to be Gaussian.
#' @param weights (New parameter for \code{hierbasis2}) Permits user-specified
#'                weights \eqn{w_{k, m}}. If left unspecified
#'                (\code{weights = NULL}) then the default is set to
#'                \eqn{w_{k, m} = k^m - (k - 1)^m}.
#'
#' @return
#' An object of class addhierbasis with the following elements
#'
#' \item{to do...}{to do...}
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
addhierbasis <- function(X, y,
                         nbasis        = 10,
                         max.lambda    = NULL,
                         lam.min.ratio = 1e-04,
                         nlam          = 50,
                         beta.mat      = NULL,
                         alpha         = NULL,
                         m.const       = 3,
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

  # each slice of the array has the orthogonal design corresponding
  # to each covariate
  design.array <- array(NA, dim = c(n, K, p))

  # initialize an empty matrix to store xbar values
  # so that we'll be able to center our data appropriately
  xbar <- matrix(NA, ncol = p, nrow = K)

  for(j in 1:p) {
    design.mat <- outer(X[, j], 1:nbasis, "^")
    if(type[1] == "gaussian") {
      design.mat.centered <- scale(design.mat, scale = FALSE)

      xbar[, j] <- attributes(design.mat.centered)[[2]]
      design.array[, , j] <- design.mat.centered
    } else {
      design.array[, , j] <- design.mat
    }

  }

  ybar <- mean(y)

  if (is.null(max.lambda)) max.lambda  <- NA
  if (is.null(alpha))      alpha       <- NA

  beta_is_zero <- all(beta.mat == 0)

  if (type[1] == "gaussian") {
    mod <- FitAdditive(y - mean(y),
                       w.mat, w,
                       design.array, beta.mat,
                       max_lambda    = max.lambda,
                       lam_min_ratio = lam.min.ratio,
                       alpha         = alpha,
                       tol           = tol,
                       p, K, n, nlam,
                       max_iter      = max.iter,
                       beta_is_zero  = beta_is_zero,
                       active_set    = colSums((beta.mat != 0) * 1))
    beta2 <- mod$beta

    # Obtain intercepts for model.
    intercept <- as.vector(ybar - (as.vector(xbar) %*% beta2))


    # Obtain the fitted values for each lambda value.
    yhats <- Matrix::crossprod(apply(design.array, 1, cbind), beta2) + ybar


    # Finally, we return an addHierbasis object.
    result <- list("beta" = beta2,
                   "intercept" = intercept,
                   "y" = y,
                   "x" = x,
                   "nbasis" = nbasis,
                   "fitted.values" = yhats,
                   "ybar" = ybar,
                   "xbar" = xbar,
                   "lam" = mod$lambdas,
                   "m.const" = m.const,
                   "type" = "gaussian")
    result$call <- match.call()

    class(result) <- "addHierBasis"

  } else {
    if(is.na(max.lambda)) {
      stop("Max lambda value needs to be specified for logistic regression.")
    }
    mod <- FitAdditiveLogistic2(y, ak.mat, ak, design.array, beta.mat,
                                intercept = intercept,
                                max_lambda = max.lambda, lam_min_ratio = lam.min.ratio,
                                alpha = alpha, tol = tol,
                                p, J, n, nlam, max_iter = max.iter,
                                beta_is_zero = beta_is_zero,
                                step_size = step.size, lineSrch_alpha = line.search.par)

    beta2 <-mod$beta

    # Obtain intercepts for model.
    intercept <- mod$intercepts


    # Obtain the fitted values for each lambda value.
    yhats <- crossprod(apply(design.array, 1, cbind), beta2)
    temp <- t(apply(yhats, 1, "+", intercept))
    phats <- 1/(1 + exp(-temp))


    # Finally, we return an addHierbasis object.
    result <- list("beta" = beta2,
                   "intercept" = intercept,
                   "y" = y,
                   "x" = x,
                   "nbasis" = nbasis,
                   "fitted.values" = phats,
                   "ybar" = ybar,
                   "xbar" = xbar,
                   "lam" = mod$lambdas,
                   "m.const" = m.const,
                   "type" = "binomial")
    result$call <- match.call()

    class(result) <- "addHierBasis"
  }

  return(result)
}
