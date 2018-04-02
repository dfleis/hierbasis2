predict.cv.additivehierbasis <- function(cv.object,
                                         new.X   = NULL,
                                         lam.idx = c("lambda.1se", "lambda.min"), ...) {
  if (is.character(lam.idx)) {
    lam.idx <- cv.object[[paste0(match.arg(lam.idx), ".idx")]]
  }

  predict(object = cv.object$model.fit, new.X = new.X, lam.idx = lam.idx, ...)
}
