#
# Computes error bars for plot.cv functions.
# Adapted from:
# https://github.com/cran/glmnet/blob/master/R/error.bars.R
#
error.bars <- function(x, upper, lower, width = 0.01, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
