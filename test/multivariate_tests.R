library(devtools)
#library(microbenchmark)
library(HierBasis)

#install_github("dfleis/hierbasis2")
library(hierbasis2)
source("R/wavelet.R")

#=========================#
#===== generate data =====#
#=========================#
set.seed(124)
n <- 500
p <- 5
sigma <- 0.5

X <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n, 0, sigma)

# define some arbitrary additive function of the covariates
# y = sum^p_{j = 1} y_j = sum^p_{j = 1} f_j(X[,j])
y1 <- X[,1] - 1.5 * X[,1]^3
y2 <- sin(pi * X[,2]) - 0.5 * sin(0.25 * pi * X[,2])
y3 <- 0.25 * X[,3]^2
for (i in 1:3) {
  y3 <- y3 + 1/i^2 * wavelet(X[,3] - 1/i^2, i)
}

ytrue <- y1 + y2 + y3
y <- ytrue + eps

plot(y1 ~ X[,1], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y2 ~ X[,2], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y3 ~ X[,3], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))

plot(y ~ X[,1], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))

#======================#
#===== fit models =====#
#======================#
pt <- proc.time()
mod1 <- HierBasis::AdditiveHierBasis(x = X, y = y, nbasis = 5)
proc.time() - pt

pt <- proc.time()
mod2 <- hierbasis2::additivehierbasis(X = X, y = y, nbasis = 5)
proc.time() - pt

yhat1 <- predict(mod1)
yhat2 <- predict(mod2)
sum(abs(yhat1 - yhat2))

pt <- proc.time()
cv.mod2 <- cv.additivehierbasis(X = X, y = y, nbasis = 10, nfolds = 10)
proc.time() - pt

m <- cv.init()

#==================#
#==== figures =====#
#==================#
ind.lam <- 40

plot(mod1, ind.func = 1, ind.lam = ind.lam, pch = 19,
     type = 'o', col = 'red', cex = 0.5, lwd = 2)
points(y1[order(X[,1])] ~ sort(X[,1]), pch = 19, cex = 0.5,
       col = rgb(0, 0, 0, 0.5))

plot(mod1, ind.func = 2, ind.lam = ind.lam, pch = 19,
     type = 'o', col = 'red', cex = 0.5, lwd = 2)
points(y2[order(X[,2])] ~ sort(X[,2]), pch = 19, cex = 0.5,
       col = rgb(0, 0, 0, 0.5))

plot(mod1, ind.func = 3, ind.lam = ind.lam, pch = 19,
     type = 'o', col = 'red', cex = 0.5, lwd = 2)
points(y3[order(X[,3])] ~ sort(X[,3]), pch = 19, cex = 0.5,
       col = rgb(0, 0, 0, 0.5))











