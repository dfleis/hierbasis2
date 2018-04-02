library(devtools)
#library(microbenchmark)

library(glmnet)
library(HierBasis)
#install_github("dfleis/hierbasis2")
library(hierbasis2)

#=====================#
#===== functions =====#
#=====================#
f1 <- function(z) 2 * z
f2 <- function(z) -1 * z^3
f3 <- function(z) 5 * sin(0.25 * pi * z)
f4 <- function(z) -2 * cos(1 * pi * z)

#=========================#
#===== generate data =====#
#=========================#
set.seed(124)
n <- 250
sigma <- 2.5

x <- sort(runif(n, -3, 3)) # sort just to make plotting easier later
eps <- rnorm(n, 0, sigma)

ytrue <- f1(x) + f2(x) + f3(x) + f4(x)
y <- ytrue + eps

plot(ytrue ~ x, ylim = range(y), pch = 19, cex = 0.25,
     col = rgb(0, 0, 0, 1), type = 'o')
points(y ~ x, pch = 19, col = rgb(0, 0, 0.75, 0.5), cex = 0.5)

#======================#
#===== fit models =====#
#======================#
pt <- proc.time()
mod.glmnet <- glmnet(x = cbind(1, x), y = y, alpha = 1, intercept = F)
proc.time() - pt

pt <- proc.time()
mod.hb1 <- HierBasis::HierBasis(x = x, y = y, nbasis = 10)
proc.time() - pt

pt <- proc.time()
mod.hb2 <- hierbasis2::hierbasis(x = x, y = y, nbasis = 10, basis.type = "wave")
proc.time() - pt

mod.hb2$fitted.values - mod.hb1$fitted.values

yhat.glmnet <- predict(mod.glmnet, newx = cbind(1, x))
yhat.hb1    <- predict(mod.hb1)
yhat.hb2    <- predict(mod.hb2)

plot(mod.hb2, sign.lambda = -1, label = F)


x.new <- rnorm(10)
yhat1 <- predict(mod.hb1, x.new)
yhat2 <- predict(mod.hb2, x.new)

sum(abs(yhat1 - yhat2))

pt <- proc.time()
cv.mod2 <- cv.hierbasis(x, y, nbasis = 10, nfolds = 10)
proc.time() - pt

new.x <- rnorm(10)
yhat1 <- predict(mod1, new.x = new.x)
yhat2 <- predict(mod2, new.x = new.x, lam.idx = cv.mod2$lambda.1se.idx)
yhat3 <- predict(cv.mod2, new.x = new.x)

yhat1[,cv.mod2$lambda.1se.idx] - yhat2
yhat2 - yhat3

cv.mod2

plot(cv.mod2)
#===================#
#===== figures =====#
#===================#

plot(ytrue ~ x, pch = 19, col = rgb(0, 0, 0, 0.5),
     cex = 0.5, ylim = range(y), type = 'o')
points(y ~ x, pch = 19, col = rgb(0, 0, 0.5, 0.5), cex = 0.75)
points(yhat2[,cv.mod2$lambda.1se.idx] ~ x, pch = 19,
       col = rgb(0, 0.75, 0, 0.35), cex = 0.5)
points(yhat2[,cv.mod2$lambda.min.idx] ~ x, pch = 19,
       col = rgb(0.75, 0, 0, 0.35), cex = 0.5)

