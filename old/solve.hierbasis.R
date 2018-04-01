# Private helper function for the univariate hierbasis
# minimization problem.
#
# Solves the penalized regression problem
#   (1/2) * || y - \Psi %*% \beta ||_2^2 + Omega(\beta),
# where
# Omega(\beta) = \sum_{k=1}^{K_n} weights[k, lam_idx] * norm(\beta[k:K_n])
# where k marks the basis function index and lam_idx tracks the
# varying values of lambda.
#
# Inputs:
#  * PSI: The n \times K centered basis expansion matrix of
#    the original covariate.
#  * y: A n-vector of centered responses.
#  * w: A n-vector of weights applying a decay to terms of
#    increasing complexity in Omega (e.g. w[k] = k^m - (k - 1)^m).
#  * w.mat: A K\times nlam matrix of weights with rows corresponding
#    to a single weight w[k] and columns corresponding to a single
#    tuning parameter lambda.
#  * max.lambda: Largest penalty.
#  * lam.min.ratio: Ratio of smallest penalty to greatest lambda.
#  * nlam: Number of penalty constants to solve the minimization
#    problem with (penalties spaced on a log scale).
#
# Outputs:
#  * betahat: A K_n \times nlam matrix of solutions to the above
#    minimization problem.
#  * lambda: The sequence of penalties used to penalize Omega.
#
.solve.hierbasis <- function(PSI, y, w, w.mat,
                             lam.min.ratio, nlam, max.lambda)
{
  n <- length(y)

  # perform scaled QR decomp. so that Q.mat^T %*% Q.mat = n * I_n
  # and Q.mat %*% R.mat = PSI
  qr.out <- qr(PSI)
  Q.mat <- qr.Q(qr.out) * sqrt(n)
  R.mat <- qr.R(qr.out) / sqrt(n)

  # See Haris et al. (2016):
  # (adapted from the HierBasis github src/hier_univariate.cpp)
  #
  # Note that for a orthogonal design with X^T %*% X = n I_n the
  # two problems are equivalent:
  #   1. (1/2n) * || Y - Q %*% beta ||^2 + lambda * Omega(beta),
  #   2. (1/2) * || Q^T %*% Y/n - beta ||^2 + lambda * Omega(beta).
  #
  # Thus all we need to do is to solve 2.

  # define v.tmp = Q^T %*% Y/n for the prox problem.
  v.tmp <- t(Q.mat) %*% y / n

  # if max.lambda value is not specified
  if (is.na(max.lambda)) {
    # TO DO... why do we use this as our upper bound?
    max.lambda <- max(abs(v.tmp) / w);
  }

  # generate lambdas
  lambda <- 10^seq(log10(max.lambda),
                   log10(max.lambda * lam.min.ratio),
                   length.out = nlam)

  # weight matrix combining w and lambda
  # to do... for some reason this flips the matrix?
  # have to transpose after apply()
  w.lam <- t(apply(w.mat, 1, function(wt) wt * lambda))

  # solve proximal problem
  betahat <- prox(v.tmp, w.lam)

  # invert of  R.mat to obtain rescale the solution back
  # to the original scale
  betahat2 <- backsolve(R.mat, betahat)

  # return results
  out        <- list()
  out$beta   <- betahat2
  out$lambda <- lambda

  return (out)
}
