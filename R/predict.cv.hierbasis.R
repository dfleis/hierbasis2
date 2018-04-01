predict.cv.hierbasis <- function(cv.object,
                                 new.x       = NULL,
                                 lam.idx     = c("lambda.1se", "lambda.min"),
                                 interpolate = FALSE, ...) {
  if (is.character(lam.idx)) {
    lam.idx <- cv.object[[paste0(match.arg(lam.idx), ".idx")]]
  }

  predict(cv.object$model.fit, new.x, lam.idx, interpolate, ...)
}
