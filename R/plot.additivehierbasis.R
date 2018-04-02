plot.additivehierbasis <- function(mod,
                                   plot.stat   = c("coef", "active"),
                                   sign.lambda = 1, ...) {

  if (plot.stat[1] == "coef") {
    plot.coef.additivehierbasis(mod, sign.lambda = sign.lambda, ...)
  } else if (plot.stat[1] == "active") {
    plot.active.additivehierbasis(mod, sign.lambda = sign.lambda, ...)
  } else {
    stop("Error in plot.additivehierbasis(): Invalid 'plot.stat'. Please enter either plot.stat = 'coef' or plot.stat = 'active'.")
  }
}
