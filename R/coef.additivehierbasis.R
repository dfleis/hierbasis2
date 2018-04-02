coef.additivehierbasis <- function(mod, lam.idx = NULL, as.arr = TRUE) {
  p      <- ncol(mod$X)
  nbasis <- mod$nbasis
  nlam   <- length(mod$lambdas)

  if (is.null(lam.idx)) {
    lam.idx <- 1:length(mod$lambdas)
  }

  beta.0 <- mod$intercept

  if (as.arr) {
    beta.X <- mod$beta.array
    attributes(beta.X)$dimnames[[1]] <- paste0("basis.", 1:nbasis)
    attributes(beta.X)$dimnames[[2]] <- paste0("X.", 1:p)
    attributes(beta.X)$dimnames[[3]] <- paste0("lam.", 1:nlam)

    names(beta.0) <- attributes(beta.X)$dimnames[[3]]

    beta.out <- list()
    beta.out$intercept <- beta.0[lam.idx]
    beta.out$X         <- beta.X[,,lam.idx]

  } else {
    beta.X <- as.matrix(mod$beta)
    beta   <- rbind(beta.0, beta.X)

    rownames(beta) <- c("Intercept",
                        paste0(rep("PSI",    times = p * nbasis), ".",
                               rep(1:p,      each = nbasis), ".",
                               rep(1:nbasis, times = p)))

    beta.out <- beta[,lam.idx]
  }

  beta.out
}
