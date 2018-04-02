print.additivehierbasis <- function(x, digits = 3, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  active.mat <- t(x$active.mat)

  out <- cbind(signif(x$lambdas, digits), x$active, active.mat)
  colnames(out) <- c("Lambda", "Active.Tot", paste0("Active.X.", 1:ncol(x$X)))
  print(out)
}
