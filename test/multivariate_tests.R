library(devtools)
#library(microbenchmark)
library(HierBasis)

#install_github("dfleis/hierbasis2")
library(hierbasis2)

#=========================#
#===== generate data =====#
#=========================#
set.seed(1)
n <- 500
p <- 6
sigma <- 1

X <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n, 0, sigma)

# define some arbitrary additive function of the covariates
# y = sum^p_{j = 1} y_j = sum^p_{j = 1} f_j(X[,j])
y1 <- 2 * sin(X[, 1])
y2 <- X[, 2]
y3 <- (X[, 3])^3/10 + sin(-0.5 * pi * X[,3])

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
mod1 <- HierBasis::AdditiveHierBasis(x = X, y = y, nbasis = 10)
proc.time() - pt

pt <- proc.time()
mod2 <- hierbasis2::additivehierbasis(X = X, y = y, nbasis = 10)
proc.time() - pt

new.X <- matrix(rnorm(5 * p), ncol = p)
yhat1 <- predict(mod1, new.X)
yhat2 <- predict(mod2, new.X)
sum(abs(yhat1 - yhat2))

pt <- proc.time()
cv.mod2 <- hierbasis2::cv.additivehierbasis(X = X, y = y, nbasis = 15, nfolds = 10)
proc.time() - pt

plot(cv.mod2)

#==================#
#==== figures =====#
#==================#
ind.lam <- cv.mod2$lambda.min.idx

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











