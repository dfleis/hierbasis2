expand.trig <- function(x, nbasis) {
  # generate and center basis exansion PSI (of order nbasis)
  # Note: The format of the trig basis expansion is changed from
  # the examples used in the original HierBasis paper.
  # 1 - Remove intercept term since we center the data.
  # 2 - Add odd degree trig functions in addition to the even
  #     degree functions.
  warning("Parameter 'basis.type = \"trig\"' not fully validated and remains experimental.")
  outer(x, 1:nbasis, FUN = function(z, nb) {
    ifelse(nb %% 2 == 1,
           cos(nb * pi * z),
           sin((nb - 1) * pi * z))
  })
}
