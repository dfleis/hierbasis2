sinc <- function(z) ifelse(z == 0, 1, sin(pi * z)/(pi * z))
wavelet <- function(z, k) {
  2^(k/2) * sinc(2^k * z)
}

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
#'               basis expansion of \code{x}. Default:
#'               \code{nbasis = length(y)}.
#' @param max.lambda The largest value of \eqn{\lambda} penalizing
#'                   the hierarchical penalty. Default:
#'                   \code{max.lambda = NULL} (computed data-adaptively).
#' @param lam.min.ratio The ratio of the smallest value of \eqn{\lambda}
#'                      to the largest \eqn{\lambda}.
#'                      Default: \code{lam.min.ratio = 1e-04}.
#' @param nlam The number \eqn{\lambda}'s to compute the \code{hierbasis}
#'             estimators (computed on a base-10 \eqn{log} scale). Default:
#'             \code{nlam = 50}.
#' @param m.const Smoothing parameter controlling the degree of polynomial
#'                smoothing/decay performed by the weights in the penalty
#'                \eqn{\Omega}.Default: \code{m.const = 3}.
#' @param type Specifies either Gaussian regression (\code{"gaussian"})
#'             with a continuous response vector or binomial regression
#'             (\code{"binomial"}) with binary response. Default type is
#'             assumed to be Gaussian.
#' @param weights (New parameter for \code{hierbasis2}) Permits user-specified
#'                weights \eqn{w_{k, m}}. If left unspecified
#'                (\code{weights = NULL}) then the default is set to
#'                \eqn{w_{k, m} = k^m - (k - 1)^m}.
#'
#' @return Returns an object of class \code{hierbasis} with elements (to finish...)
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
                      m.const       = 3,
                      type          = c("gaussian", "binomial"),
                      weights       = NULL,
                      basis.type    = c("poly", "trig", "wave"))
{
  # To do:
  #   * Implement binomial hierbasis estimator.
  #   * Consider multicategory responses (not in HierBasis).
  #   * Add max.iter, tol, max.iter.inner, tol.inner parameters (only
  #     used in the 'binomial' case).
  #   * Implementing nonpolynomial basis expansion (wavelet, trig, etc).
  #   * Implementing other weight structures.

  # extract number of observations
  n <- length(y)

  if (basis.type[1] == "poly") {
    # POLY basis expansion

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
    # TRIG basis expansion

    # generate and center basis exansion PSI (of order nbasis)
    PSI <- outer(x, 1:nbasis, FUN = function(z, nb) {
      ifelse(nb == 1, rep(1, length(z)), {
        kb <- floor(nb/2) * 2
        ifelse(nb %% 2 == 0,
               cos(kb * pi * z),
               sin(kb * pi * z))
      })
    })
    PSI.c <- scale(PSI, scale = F)
    PSIbar <- attributes(PSI.c)[[2]]

    warning("Warning in hierbasis2::hierbasis(). Parameter 'basis.type = \"trig\"'
             is unfinished and experimental.")

  } else if (basis.type[1] == "wave") {
    # WAVE basis expansion

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

  # center responses
  ybar <- mean(y)
  y.c <- y - ybar

  # compute penalty weights
  if (is.null(weights)) {
    w <- (1:nbasis)^m.const - (0:(nbasis - 1))^m.const
  } else {
    # confirm weights has the correct length of nbasis
    if (length(weights) != nbasis) {
      warning("Warning in hierbasis2::hierbasis(): 'weights' vector should be length 'nbasis'.")
    }
    w <- weights
  }

  # define weight matrix data structrue (which will include
  # the tuning parameter lambda)
  # rows correspond to a single weight w_k value
  # columns corresponds to a single tuning parameter lambda
  w.mat <- matrix(rep(w, nlam), ncol = nlam)

  if (is.null(max.lambda)) {
    max.lambda <- NA
  }

  # create empty the hierbasis object to be returned
  out <- list(); class(out) <- "hierbasis"

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
    out$m.const               <- m.const
    out$nbasis                <- nbasis
    out$basis.expansion       <- PSI
    out$basis.expansion.means <- PSIbar
    out$ybar                  <- ybar

    out$lambdas       <- as.vector(hierbasis.fit$lambda)
    out$intercept     <- intercept
    out$beta          <- betahat
    out$fitted.values <- t(yhat)
    out$active        <- active.set
    out$basis.type    <- basis.type
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

print.hierbasis <- function(x, digits = 3, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Lambda = signif(x$lambdas, digits),
              Deg.of.Poly = x$active))
}

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
                              interpolate = FALSE, ...) {
  nlam <- length(object$lambdas)  # Number of lambdas.

  if (!interpolate) {
    if (is.null(new.x)) {
      # return the same fitted response if no new predictors
      # are provided
      object$fitted.values
    } else {
      if (object$basis.type == "poly") {
        # POLY: basis-expanded design matrix

        newx.mat <- sapply(1:object$nbasis, function(i) {
          new.x^i
        })
      } else if (object$basis.type == "trig") {
        # TRIG: basis-expanded design matrix

        newx.mat <- outer(x, 1:object$nbasis, FUN = function(z, nb) {
          ifelse(nb == 1, rep(1, length(z)), {
            kb <- floor(nb/2) * 2
            ifelse(nb %% 2 == 0,
                   cos(kb * pi * z),
                   sin(kb * pi * z))
          })
        })

      } else if (object$basis.type == "wave") {
        # WAVE: basis-expanded design matrix

        newx.mat <- outer(x, 0:(object$nbasis - 1), FUN = function(z, nb) {
          wavelet(z, nb)
        })

      }

      # X %*% beta without the intercept
      fitted <- newx.mat %*% object$beta
      # add the intercept
      t(apply(fitted, 1, "+", object$intercept))
    }
  } else {
    if (is.null(new.x)) {
      return (object$fitted.values)
    }

    # Return predicted values.
    sapply(1:nlam, FUN = function(i) {
      # Obtain curve for a particular value.
      yhat.temp <- object$fitted.values[, i]

      # Return predictions.
      approx(x = object$x, y = yhat.temp, xout = new.x)$y
    })

  }
}









