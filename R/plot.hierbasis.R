plot.hierbasis <- function(mod,
                           label       = FALSE,
                           sign.lambda = 1, ...) {
  nlam   <- length(mod$lambdas)
  nbasis <- mod$nbasis

  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))

  plot.args <- list()
  plot.args$x    <- sign.lambda * log(mod$lambdas)
  plot.args$xlab <- xlab
  plot.args$type <- "n"
  plot.args$y    <- mod$beta
  plot.args$ylab <- "Coefficients"

  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args

  top.x <- pretty(plot.args$x, shrink.sml = 5)
  top.labs <- approx(x = plot.args$x, y = mod$active, xout = top.x,
                     rule = 2, method = "constant", f = 1)

  xlim <- range(plot.args$x)
  ylim <- range(plot.args$y)
  if (label) {
    ylim.center <- mean(ylim)
    ylim.diff   <- diff(ylim)
    if (sign.lambda > 0) {
      xlim[1] <- xlim[2] - diff(xlim) * 1.1
    } else {
      xlim[2] <- xlim[1] + diff(xlim) * 1.1
    }
    ylim[1] <- ylim.center - ylim.diff/2*1.05
    ylim[2] <- ylim.center + ylim.diff/2*1.05
  }

  plot(NA,
       xlim = xlim,
       ylim = ylim,
       xlab = plot.args$xlab,
       ylab = plot.args$ylab)
  matplot(plot.args$x, t(plot.args$y), type = 'l',
          pch = NULL,
          lwd = 1.5,
          add = T)
  axis(3, at = top.labs$x, labels = top.labs$y, tcl = NA)


  if (label) {
    pos <- ifelse(sign.lambda > 0, 2, 4)

    for (k in 1:nrow(mod$beta)) {
      text(x = plot.args$x[nlam], y = plot.args$y[k, nlam],
           labels = substitute(paste(beta[k]), list(k = k)),
           cex = 0.75, pos = pos)
    }
  }
}

