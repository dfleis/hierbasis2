#
# Adapted from:
# https://github.com/cran/glmnet/blob/master/R/plot.cv.glmnet.R
#
plot.cv.hierbasis <- function(cv.obj, sign.lambda = 1, ...) {
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))

  plot.args <- list(
    x = sign.lambda * log(cv.obj$lambdas),
    y = cv.obj$test.err,
    ylim = range(cv.obj$test.err.lo, cv.obj$test.err.hi),
    xlab = xlab,
    ylab = cv.obj$loss.func,
    type = "n"
  )
  new.args <- list(...)

  if (length(new.args)) plot.args[names(new.args)] <- new.args

  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cv.obj$lambdas),
             cv.obj$test.err.hi,
             cv.obj$test.err.lo,
             col = 'gray60')
  points(sign.lambda * log(cv.obj$lambdas), cv.obj$test.err,
         pch = 20, col = 'red')
  axis(side   = 3,
       at     = sign.lambda * log(cv.obj$lambdas),
       labels = paste(cv.obj$model.fit$active),
       tick   = F,
       line   = 0)
  abline(v = sign.lambda * log(cv.obj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cv.obj$lambda.1se), lty = 3)
  invisible()
}
