expand.wave <- function(x, nbasis) {
  # wavelet expansion
  warning("Parameter 'basis.type = \"wave\"' not fully validated and remains experimental.")
  PSI <- outer(x, 0:(nbasis - 1), FUN = function(z, nb) {
    wavelet(z, nb)
  })
}
