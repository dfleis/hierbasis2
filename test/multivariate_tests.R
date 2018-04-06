library(devtools)
#library(microbenchmark)

library(glmnet)
library(HierBasis)
#install_github("dfleis/hierbasis2")
library(hierbasis2)


#=========================#
#===== generate data =====#
#=========================#
set.seed(1)
n <- 250
p <- 3
sigma <- 0

X <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n, 0, sigma)

# define some arbitrary additive function of the covariates
# y = sum^p_{j = 1} y_j = sum^p_{j = 1} f_j(X[,j])
y1 <- 2 * sin(X[, 1])
y2 <- X[, 2]
y3 <- 0.25 * X[,3]^2 - 3 * exp(-X[,3]^10) + cos(X[,3])

ytrue <- y1 + y2 + y3
y <- ytrue + eps

plot(y1 ~ X[,1], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y2 ~ X[,2], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y3 ~ X[,3], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))

plot(y ~ X[,3], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))

#======================#
#===== fit models =====#
#======================#
pt <- proc.time()
mod.glmnet <- glmnet(x = X, y = y, alpha = 1)
proc.time() - pt

pt <- proc.time()
mod.ahb1 <- HierBasis::AdditiveHierBasis(x = X, y = y, nbasis = 20)
proc.time() - pt

pt <- proc.time()
mod.ahb2 <- hierbasis2::additivehierbasis(X = X, y = y, nbasis = 20)
proc.time() - pt

pt <- proc.time()
cv.mod2 <- hierbasis2::cv.additivehierbasis(X = X, y = y, nbasis = 20, nfolds = 10)
proc.time() - pt

new.X <- matrix(rnorm(9 * p), nrow = 9)
yh1 <- predict(mod.ahb1, new.x = new.X)
yh2.1 <- predict(cv.mod2, new.X = new.X, lam.idx = "lambda.1se")
yh2.2 <- predict(cv.mod2, new.X = new.X, lam.idx = "lambda.min")

plot(cv.mod2)
plot(cv.mod2$test.err[1:20])
points(cv.mod2$test.err.hi[1:20], col = 'red')

yh2.1 - yh2.2

print(mod.ahb2)
plot(mod.ahb2,
     pred.idx    = 2,
     sign.lambda = -1,
     plot.stat   = "active",
     plot.type   = "image",
     legend      = T)

plot(mod.ahb2,
     plot.stat   = "active",
     plot.type   = "image",
     legend      = T)


yhat <- predict(cv.mod2)
str(yhat)


y
yh


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











