plot.active.additivehierbasis <- function(mod,
                                          sign.lambda = 1,
                                          plot.type   = c("image", "lines"),
                                          col.image   = parula(64),
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
  plot.args$main <- "Active Set"

  if (plot.type[1] == "image") {
    active.flip <- t(plot.args$y)
    xax         <- plot.args$x
    if (sign.lambda > 0) {
      active.flip <- apply(active.flip, 2, rev)
      xax         <- rev(xax)
    }
    image.mat(x      = xax,
              y      = 1:p,
              z      = active.flip,
              legend = legend,
              col    = col.image,
              xlab   = plot.args$xlab,
              ylab   = "Covariate",
              main   = plot.args$main)

  } else if (plot.type[1] == "lines") {
    top.x <- pretty(plot.args$x, n = length(pretty(plot.args$x)) + 2)
    top.labs <- approx(x = plot.args$x, y = mod$active,
                       xout   = top.x,
                       rule   = 3,
                       method = "constant",
                       f      = 1)
    plot(NA,
         xlim = range(plot.args$x),
         ylim = range(plot.args$y),
         xlab = plot.args$xlab,
         ylab = plot.args$ylab,
         main = plot.args$main)
    matplot(plot.args$x, t(plot.args$y), type = 'l',
            pch = NULL,
            lwd = 1.5,
            add = T)
    axis(3, at = top.labs$x, labels = top.labs$y, tcl = NA, padj = 1)

  } else {
    stop("Error in plot.active.additivehierbasis(): Invalid 'plot.type' specified. Only defined for 'image' and 'lines'.")
  }
}
