plot.active.additivehierbasis <- function(mod,
                                          sign.lambda = 1,
                                          plot.type   = c("image", "lines"),
                                          legend      = F, ...) {
  p          <- ncol(mod$X)
  lambdas    <- mod$lambdas
  active.mat <- mod$active.mat

  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))

  plot.args <- list()
  plot.args$x    <- sign.lambda * log(lambdas)
  plot.args$xlab <- xlab
  plot.args$y    <- active.mat
  plot.args$ylab <- "Active"

  if (plot.type[1] == "image") {
    active.flip <- t(plot.args$y)
    xax         <- plot.args$x
    if (sign.lambda > 0) {
      active.flip <- apply(active.flip, 2, rev)
      xax         <- rev(xax)
    }

    if (legend) {
      image.plot(x    = xax,
                 y    = 1:p,
                 z    = active.flip,
                 col  = parula(64),
                 xlab = plot.args$xlab,
                 ylab = "Covariate")
    } else {
      image(x    = xax,
            y    = 1:p,
            z    = active.flip,
            col  = parula(64),
            xlab = plot.args$xlab,
            ylab = "Covariate")
    }

  } else if (plot.type[1] == "lines") {
    plot(NA,
         xlim = range(plot.args$x),
         ylim = range(plot.args$y),
         xlab = plot.args$xlab,
         ylab = plot.args$ylab)
    matplot(plot.args$x, t(plot.args$y), type = 'l',
            pch = NULL,
            lwd = 1.5,
            add = T)
  } else {
    stop("Error in plot.active.additivehierbasis(): Invalid 'plot.type' specified. Only defined for 'image' and 'lines'.")
  }



}
