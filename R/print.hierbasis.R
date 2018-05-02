print.hierbasis <- function(x, digits = 3, ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Lambda = signif(x$lambdas, digits),
              Deg.of.Poly = x$active))
}
