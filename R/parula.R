#' Parula color map
#'
#' A useful color map for image plots.
#' @param n Number of color levels. Default setting \code{n = 64}.
#' @export
parula <- function(n = 64) {
  colormap <- colorRampPalette(
  c(
    rgb(53,  42,  135, maxColorValue = 255),
    rgb(15,  92,  221, maxColorValue = 255),
    rgb(18,  125, 216, maxColorValue = 255),
    rgb(7,   156, 207, maxColorValue = 255),
    rgb(21,  177, 180, maxColorValue = 255),
    rgb(89,  189, 140, maxColorValue = 255),
    rgb(165, 190, 107, maxColorValue = 255),
    rgb(225, 185, 82,  maxColorValue = 255),
    rgb(252, 206, 46,  maxColorValue = 255),
    rgb(249, 251, 14,  maxColorValue = 255)
  ))
  return (colormap(n))
}
