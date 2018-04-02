plot.coef.additivehierbasis <- function(mod,
                                        pred.idx,
                                        sign.lambda = 1,
                                        plot.type   = c("image", "lines"),
                                        col.image   = redblue.grad(65),
                                        legend      = F, ...) {
  # reblue.grad() should be given an odd number so that the
  # midpoint will be associated with the white midpoint
  # between 'red' and 'blue
  if (length(pred.idx) != 1) {
    stop("Error in plot.coef.additivehierbasis(): Parameter 'pred.idx' must be a single positive integer in the range of 0, 1, ..., p, where pred.idx = 0 corresponds to the intercept.")
  }
  if (pred.idx == 0 & length(plot.type) == 2) {
    plot.type <- "lines"
  }

  p       <- ncol(mod$X)
  nbasis  <- mod$nbasis
  lambdas <- mod$lambdas

  beta.arr <- coef(mod)
  if (pred.idx == 0) {
    beta <- beta.arr$intercept
  } else {
    beta   <- beta.arr$X[,pred.idx,]
    active <- mod$active.mat[pred.idx,]
  }

  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))

  plot.args <- list()
  plot.args$main <- ifelse(pred.idx == 0, "Intercept", paste0("Covariate ", pred.idx))
  plot.args$x    <- sign.lambda * log(lambdas)
  plot.args$xlab <- xlab
  plot.args$y    <- beta
  plot.args$ylab <- expression(beta)

  if (plot.type[1] == "image") {
    if (pred.idx == 0) {
      stop("Error in plot.coef.additivehierbasis(): plot.type = 'image' not defined for intercept estimates. Please select plot.type = 'lines' or select pred.idx to a positive integer in 1, ..., p.")
    }

    plot.args$main <- paste0(plot.args$main, ": ", "Coefficients")

    beta.flip <- t(plot.args$y)
    xax       <- plot.args$x
    if (sign.lambda > 0) {
      beta.flip <- apply(beta.flip, 2, rev)
      xax       <- rev(xax)
    }

    zm <- max(abs(beta.flip))
    zl <- c(-zm, zm)

    image.mat(x      = xax,
              y      = 1:nbasis,
              z      = beta.flip,
              legend = legend,
              col    = col.image,
              xlab   = plot.args$xlab,
              ylab   = "Basis Degree",
              zlim   = zl,
              main   = plot.args$main)

  } else if (plot.type[1] == "lines") {

    plot(NA,
         xlim = range(plot.args$x),
         ylim = range(plot.args$y),
         xlab = plot.args$xlab,
         ylab = plot.args$ylab,
         main = plot.args$main)

    if (pred.idx > 0) {
      top.x <- pretty(plot.args$x, n = length(pretty(plot.args$x)) + 2)
      top.labs <- approx(x = plot.args$x, y = active,
                         xout   = top.x,
                         rule   = 3,
                         method = "constant",
                         f      = 1)
      matplot(plot.args$x, t(plot.args$y), type = 'l',
              pch = NULL,
              lwd = 1.5,
              add = T)
      axis(3, at = top.labs$x, labels = top.labs$y, tcl = NA, padj = 1)
    } else {
      lines(plot.args$x, plot.args$y, lwd = 1.5)
    }
  } else {
    stop("Error in plot.coef.additivehierbasis(): Invalid 'plot.type' specified. Only defined for 'image' and 'lines'.")
  }
}








