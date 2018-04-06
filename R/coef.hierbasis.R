coef.hierbasis <- function(mod, lam.idx = NULL) {
  if (is.null(lam.idx)) {
    lam.idx <- 1:length(mod$lambdas)
  }
  coef.out <- as.matrix(rbind(mod$intercept, mod$beta)[,lam.idx])
  rownames(coef.out) <- c("Intercept", paste0("PSI.", 1:mod$nbasis))
  colnames(coef.out) <- paste0("lam.", lam.idx)
  coef.out
}
