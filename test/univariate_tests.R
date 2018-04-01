library(devtools)
#library(microbenchmark)
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
sigma <- 0.5

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
mod1 <- HierBasis::HierBasis(x = x, y = y, nbasis = 10)
proc.time() - pt

pt <- proc.time()
mod2 <- hierbasis2::hierbasis(x = x, y = y, nbasis = 10)
proc.time() - pt


mod2$nbasis
x.new <- rnorm(10)
yhat1 <- predict(mod1, x.new)
yhat2 <- predict(mod2, x.new)

sum(abs(yhat1 - yhat2))



pt <- proc.time()
cv.mod2 <- cv.hierbasis(x, y, nbasis = 10, nfolds = 10)
proc.time() - pt

yhat1 <- predict(mod1)
yhat2 <- predict(mod2)
yhat2 <- predict(cv.mod2$model.fit)

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

