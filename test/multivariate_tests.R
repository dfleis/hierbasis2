library(devtools)
#library(microbenchmark)
library(HierBasis)

#install_github("dfleis/hierbasis2")
library(hierbasis2)

#=========================#
#===== generate data =====#
#=========================#
set.seed(124)
n <- 250
p <- 5
sigma <- 0.5


X <- matrix(rnorm(n * p), nrow = n)
beta <- 1:p
eps <- rnorm(n, 0, sigma)
y <- X %*% beta + eps

#======================#
#===== fit models =====#
#======================#
pt <- proc.time()
mod1 <- HierBasis::AdditiveHierBasis(x = X, y = y, nbasis = 10)
proc.time() - pt





