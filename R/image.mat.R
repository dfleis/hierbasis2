image.mat <- function(x, y, z, legend, ...) {
  if (legend) {
    image.plot(x = x, y = y, z = z, ...)
  } else {
    image(x = x, y = y, z = z, ...)
  }
}

