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
n <- 500
p <- 6
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

print(mod.ahb2)
plot.active.additivehierbasis(mod.ahb2, sign.lambda = -1)

new.X <- matrix(rnorm(10 * p), ncol = p)
yhat1 <- predict(mod1, new.X)
yhat2 <- predict(mod2, new.X = new.X, lam.idx = cv.mod2$lambda.1se.idx)
yhat1[,cv.mod2$lambda.1se.idx] - yhat2

yh2.0 <- predict(mod2, new.X = new.X, lam.idx = cv.mod2$lambda.1se.idx)
yh2.1 <- predict(cv.mod2, new.X = new.X, lam.idx = cv.mod2$lambda.1se.idx)
yh2.2 <- predict(cv.mod2, new.X = new.X)
yh2.0 - yh2.2

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











