#
# Helper function for cv.hierbasis
# (and possibly cv.addhierbasis, but I havent gotten to it yet).
#
# Finds the lambdas associated with the minimum and 1-standard-error
# cross-validation testing error.
#
# Adapted from https://github.com/cran/glmnet/blob/master/R/getmin.R
#
getmin.lambda <- function(lambdas, cv.mean, cv.se) {
  cv.min      <- min(cv.mean, na.rm = T)
  cv.mean.idx <- cv.mean <= cv.min

  # Note: Instead of taking the max over all lambdas satisfying
  # cv.mean <= cv.min (i.e., take the first lambda corresponding
  # to the indices satisfying cv.mean <= cv.min) we COULD have just
  # done which(cv.mean == cv.min). I'm not 100% sure why it was done
  # like this, but I think we can understand it in a few ways:
  # 1 - We need a way of choosing a lambda in the unlikely event of
  #     2+ cv.mean errors being tied. Choosing the max lambda
  #     satisfying this inequality implies a more conservative estimate
  #     of our coefficients and fitted responses (a larger lambda
  #     leads to a solution closer to the mean (ybar) and will be less
  #     likely to lead to overfitting compared to smaller lambdas).
  # 2 - Maybe this method is faster than explicitly selecting the
  #     minimizing lambda (calling which(cv.mean == cv.min) may not
  #     be as speedy as computing the max).
  lambda.min <- max(lambdas[cv.mean.idx], na.rm = T)
  lambda.min.idx <- match(lambda.min, lambdas)

  cv.1se.min <- (cv.mean + cv.se)[lambda.min.idx]
  cv.1se.idx <- cv.mean <= cv.1se.min
  lambda.1se <- max(lambdas[cv.1se.idx], na.rm = T)
  lambda.1se.idx <- match(lambda.1se, lambdas)

  list(
    "lambda.min"     = lambda.min,
    "lambda.min.idx" = lambda.min.idx,
    "lambda.1se"     = lambda.1se,
    "lambda.1se.idx" = lambda.1se.idx)
}


