coef.cv.additivehierbasis <- function(cv.mod,
                                      lam.idx = c("lambda.1se", "lambda.min"), ...) {
  if (is.character(lam.idx)) {
    lam.idx <- cv.mod[[paste0(match.arg(lam.idx), ".idx")]]
  }
  coef(cv.mod$model.fit, lam.idx, ...)
}

