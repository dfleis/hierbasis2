nu <- function(x) {
  ifelse(x <= 0, 0,
         ifelse(x >= 1, 1, x))
}
wavelet.meyer.scale <- function(omega) {
  1/sqrt(2 * pi) *
    ifelse(abs(omega) < 2 * pi / 3, 1,
           ifelse(abs(omega) > 4 * pi / 3, 0,
                  cos(pi/2 * nu(3 / (2 * pi) * abs(omega) - 1))
                  )
           )
}
wavelet.meyer.spec <- function(omega) {
  I <- complex(real = 0, imaginary = 1)
  1/sqrt(2 * pi) *
    ifelse(abs(omega) < 2 * pi/3, 0,
           ifelse(abs(omega) <= 4 * pi/3,
                  sin(pi/2 * nu(3 * abs(omega)/(2 * pi) - 1)) * exp(I * omega/2),
                  ifelse(abs(omega) <= 8 * pi/3,
                         cos(pi/2 * nu(3 * abs(omega)/(4 * pi) - 1)) * exp(I * omega/2),
                         0)
           )
    )
}


z <- seq(-10, 10, length.out = 1e3)
plot(wavelet.meyer.scale(z) ~ z, type = 'l', lwd = 2, col = 'darkblue')
plot(abs(wavelet.meyer.spec(z)) ~ z, type = 'l', lwd = 2, col = 'darkblue')

plot(wavelet(z) ~ z, type = 'l', lwd = 2, col = 'darkblue')

z <- seq(0, 1, length.out = 1e3)
plot(cos(2 * pi * z) ~ z, type = 'l', lwd = 2, col = 'darkblue')
plot(sin(2 * pi * z) ~ z, type = 'l', lwd = 2, col = 'darkblue')
plot(cos(4 * pi * z) ~ z, type = 'l', lwd = 2, col = 'darkblue')
