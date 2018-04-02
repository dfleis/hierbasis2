library(devtools)
#library(microbenchmark)
#install_bitbucket("researchAH/HierBasisRpackage")
install_github("asadharis/HierBasis")
library(HierBasis)

#install_github("dfleis/hierbasis2")
library(hierbasis2)

#=====================#
#===== functions =====#
#=====================#
logit <- function(p) log(p/(1 - p))
expit <- function(z) 1/(1 + exp(-z))

#=========================#
#===== generate data =====#
#=========================#
set.seed(124)
n <- 250

x <- sort(rnorm(n))
eta <- x
pi <- expit(eta)
y  <- rbinom(n = n, size = 1, prob = pi)

mod.glm <- glm(y ~ x, family = binomial)

etahat <- predict(mod.glm)
plot(y ~ x)
lines(expit(etahat) ~ x, type = 'l', col = 'red')



#======================#
#===== fit models =====#
#======================#
pt <- proc.time()
mod1 <- HierBasis::HierBasis(x = x, y = y,
                             nbasis         = 20,
                             nlam           = 10,
                             type           = "binomial",
                             max.iter       = 100,
                             tol            = 1e-3,
                             max.iter.inner = 100,
                             tol.inner      = 10)
proc.time() - pt

plot(y ~ x, ylim = range(y), pch = 19, cex = 0.5,
     col = rgb(0, 0, 0, 1))
for (k in 1:length(mod1$lambdas)) {
  points(mod1$fitted.values[,k] ~ x, col = 'red', cex = 0.25)
}









