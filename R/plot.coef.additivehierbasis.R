plot.coef.additivehierbasis <- function(mod,
                                        sign.lambda = 1,
                                        legend      = F,
                                        ...) {
  # beta.arr  <- coef(mod)
  # beta.0    <- beta$intercept
  # beta.pred <- beta$X

  beta <- coef(mod, as.arr = F)
  str(beta)

  image.mat(x = 1:121, y = 1:50, z = beta, legend = F)
}
