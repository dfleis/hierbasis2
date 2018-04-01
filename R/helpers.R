.sinc <- function(z) ifelse(z == 0, 1, sin(pi * z)/(pi * z))

.wavelet <- function(z, k) {
  2^(k/2) * .sinc(2^k * z)
}
