#' Cross-validation for Multivariate additivehierbasis Models
#'
#' @details The function runs \code{nfolds} + 1 times. If \code{lambdas}
#'          is not specified then it uses the first run to find the
#'          sequence of \code{lambdas}. The remaining \code{nfolds}
#'          runs are done using the fold subsets as usual.
#'
#' @param X An \eqn{n \times p}{n x p} matrix of covariates.
#' @param y A univariate vector respresenting the response.
#' @param lambda (Optional) User-specified sequence of tuning parameters \eqn{lambda}.
#' @param nfolds Number of cross-validation folds. Default: \code{nfolds = 10}.
#' @param ... Other arguments that may be passed to \code{hierbasis}.
#'
#' @return Returns an object of class \code{hierbasis} with elements (to finish...)
#'
#' \item{to do...}{to do...}
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
cv.additivehierbasis <- function(X, y, lambdas = NULL, nfolds = 10, ...) {
  # TO DO: Allow for other loss functions than rmse.
  # TO DO: Allow for parallel computing (foreach library?).
  # TO DO: Split some overhead off in its own function.

  n <- length(y)
  mod.init <- cv.init(x = X, y = y, lambdas = lambdas, ...)

  # assign observations to cross-validation
  # fold 1, ..., nfolds
  folds <- sample(cut(1:n, breaks = nfolds, labels = F))

  # get training subset indices
  train.idx <- lapply(1:nfolds, function(i) !(folds %in% i))

  cv.err <- sapply(train.idx, function(trn) {
    tst <- !trn # test subset indices

    X.trn <- X[trn,]; y.trn <- y[trn]
    X.tst <- X[tst,]; y.tst <- y[tst]

    # compute training model and training error
    mod.cv <- additivehierbasis(X = X.trn, y = y.trn,
                                nlam          = mod.init$nlam,
                                max.lambda    = mod.init$max.lambda,
                                lam.min.ratio = mod.init$lam.min.ratio, ...)

    yhat.trn <- mod.cv$fitted.values # extract training fits
    yhat.tst <- predict(mod.cv, new.X = X.tst) # compute test fits

    # compute errors
    if (mod.cv$type[1] == "binomial") {
      warning("Warning: 'type' not yet defined for 'binomial' models.")
    } else {
      trn.err <- apply(yhat.trn, 2, function(yhat) rmse(y.trn, yhat))
      tst.err <- apply(yhat.tst, 2, function(yhat) rmse(y.tst, yhat))
    }

    list("train" = trn.err, "test" = tst.err)
  })

  cv.train.err <- do.call("cbind", cv.err["train",])
  cv.test.err  <- do.call("cbind", cv.err["test",])

  # average over all folds
  # each entry is an error corresponding to a unique lambda
  train.err <- rowSums(cv.train.err)/nfolds
  test.err  <- rowSums(cv.test.err)/nfolds
  # compute standard errors
  test.err.se <- apply(cv.test.err, 1, function(z) sd(z)/sqrt(length(z)))
  test.err.hi <- test.err + test.err.se
  test.err.lo <- test.err - test.err.se

  best.lambdas <- getmin.lambda(mod.init$lambdas, test.err, test.err.se)

  out <- list()
  out$additivehierbasis.fit <- mod.init$model.fit$mod.init
  out$train.err             <- train.err
  out$test.err              <- test.err
  out$test.err.se           <- test.err.se
  out$test.err.hi           <- test.err.hi
  out$test.err.lo           <- test.err.lo
  out$lambdas               <- mod.init$lambdas
  out$loss.func             <- "rmse" # to do... generalize...
  out <- c(out, best.lambdas)

  class(out) <- "cv.additivehierbasis"
  out
}
