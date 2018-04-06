library(foreach)
nfolds <- 1e3
parallel <- TRUE

n <- 1e3
pt <- proc.time()
sim1 <- foreach(i = seq(nfolds), .packages = c('hierbasis2')) %dopar% {
  rnorm(n)
}
proc.time() - pt

pt <- proc.time()
sim2 <- replicate(nfolds, rnorm(n))
proc.time() - pt

