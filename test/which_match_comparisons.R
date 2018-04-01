#
# Compares the speed of different techniques finding the
# elements in a vector
#
# 1 - which(x == x_i)
# 2 - match(x_i, x)
# 3 - sum((1:length(x)) * (x == 0))

library(microbenchmark)

#===== functions =====#
my_check <- function(vals) {
  # function for use in microbenchmark()
  # to confirm all functions are indeed returning
  # the same output
  all(
    sapply(vals[-1],
         function(v) {
           identical(vals[[1]], v)
           }
         )
  )
}

f1 <- function(x, x_i) {
  which(x == x_i)
}
f2 <- function(x, x_i) {
  match(x_i, x)
}

FUNS <- list("f1" = f1, "f2" = f2)

#===== parameters =====#
n.max    <- log10(1e4)
n.min    <- log10(1e3)
n.len    <- 10
times.mb <- 50
unit.mb  <- "us"

#===== generate data =====#
n.samp <- floor(10^seq(n.min, n.max, length.out = n.len))
x.list <- sapply(n.samp, FUN = function(n.s) sample(c(0, rnorm(n.s - 1))))

time.out.colnames <- c("min", "lq", "mean", "median", "uq", "max")
time.out <- array(NA, dim = c(length(FUNS),
                              length(time.out.colnames),
                              length(n.samp)))
attributes(time.out)$dimnames[[1]] <- names(FUNS)
attributes(time.out)$dimnames[[2]] <- time.out.colnames

time.ratios <- vector(mode = 'list', length = length(n.samp))

#==== do work: time code =====#

for (i in 1:length(n.samp)) {
  x <- sample(x.list[[i]])
  pt <- proc.time()
  mb <- microbenchmark(f1(x, 0), f2(x, 0),
                       times = times.mb, unit = unit.mb, check = my_check)
  tm.tmp <- proc.time() - pt

  smb <- summary(mb)
  time.out[,,i] <- as.matrix(smb[,colnames(time.out)])

  time.ratios[[i]] <- mb$time[mb$expr == levels(mb$expr)[1]]/mb$time[mb$expr == levels(mb$expr)[2]]
  print(c(i, tm.tmp))
}

#===== figures =====#
fun_cols <- rainbow(length(FUNS))
fun_ltys <- c("dashed", "solid", "dotdash")

plot(NA, log = 'xy',
     xlim = range(n.samp),
     ylim = range(time.out[, "min",], 1.25 * time.out[, "uq",]),
     xlab = "Nb. of Samples",
     ylab = substitute(paste("Time (", mu, "s)")))
for (j in 1:length(FUNS)) {
  segments(x0 = n.samp, x1 = n.samp, lwd = 1.25,
           y0 = time.out[j, "lq",],
           y1 = time.out[j, "uq",],
           col = fun_cols[j])
}
for (j in 1:length(FUNS)) {
  points(time.out[j, "median",] ~ n.samp, type = 'o',
         pch = 21, bg = 'white', cex = 0.75, lwd = 2,
         col = fun_cols[j],
         lty = fun_ltys[j])
}
legend("topleft", legend = names(FUNS), bty = 'n',
       col = fun_cols, lty = fun_ltys, lwd = 1.5)


time.ratios.mean <- sapply(time.ratios, mean)
time.ratios.q <- sapply(time.ratios, quantile, probs = c(0.25, 0.5, 0.75))

plot(NA, log = 'xy',
     xlim = range(n.samp),
     ylim = c(0.5 * min(time.ratios.q),
              1.25 * max(time.ratios.q)),
     xlab = "Nb. of Samples",
     ylab = "(which() time)/(match() time)")
abline(h = 1, lwd = 0.75, lty = 'dashed', col = 'gray60')
abline(v = 16252, col = 'red', lwd = 2)
segments(x0 = n.samp, x1 = n.samp,
         y0 = time.ratios.q["75%",],
         y1 = time.ratios.q["25%",],
         col = 'darkblue', lwd = 1)
points(time.ratios.q["50%",] ~ n.samp, type = 'o',
       pch = 21, bg = 'white', lwd = 1.5, col = 'darkblue')




