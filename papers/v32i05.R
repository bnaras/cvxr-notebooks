###################################################
### chunk number 1: preliminaries
###################################################
library("isotone")


###################################################
### chunk number 2: 
###################################################
require("isotone")
data("pituitary")
head(pituitary)


###################################################
### chunk number 3: 
###################################################
res1 <- with(pituitary, gpava(age, size, ties = "primary"))
res2 <- with(pituitary, gpava(age, size, ties = "secondary"))
res3 <- with(pituitary, gpava(age, size, ties = "tertiary"))


###################################################
### chunk number 4: pit-plot eval=FALSE
###################################################
## layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE))
## plot(res1, main = "PAVA plot (primary)")
## plot(res2, main = "PAVA plot (secondary)")
## plot(res3,  main = "PAVA plot (tertiary)")


###################################################
### chunk number 5: pit-plot1
###################################################
layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = TRUE))
plot(res1, main = "PAVA plot (primary)")
plot(res2, main = "PAVA plot (secondary)")
plot(res3,  main = "PAVA plot (tertiary)")


###################################################
### chunk number 6: 
###################################################
tapply(res3$x, res3$z, mean)


###################################################
### chunk number 7: 
###################################################
data("posturo")
head(posturo)


###################################################
### chunk number 8: 
###################################################
res.mean <- with(posturo, gpava(height, cbind(SOT.1, SOT.2, SOT.3),
  solver = weighted.mean, ties = "secondary"))
res.median <- with(posturo, gpava(height, cbind(SOT.1, SOT.2, SOT.3),
  solver = weighted.median, ties = "secondary"))


###################################################
### chunk number 9: post-plot eval=FALSE
###################################################
## plot(res.mean)
## plot(res.median)


###################################################
### chunk number 10: post-plot1
###################################################
par(mar = c(5,4,4,2))
par(mfrow = c(1, 2))
plot(res.mean)
plot(res.median)


###################################################
### chunk number 11: 
###################################################
set.seed(12345)
y <- rnorm(9)
w1 <- rep(1, 9)
Atot <- cbind(1:8, 2:9)


###################################################
### chunk number 12: 
###################################################
Atot


###################################################
### chunk number 13: 
###################################################
fit.ls1 <- activeSet(Atot, "LS", y = y, weights = w1)
fit.ls2 <- activeSet(Atot, fSolver, y = y, weights = w1,
  fobj = function(x) sum(w1 * (x - y)^2),
  gobj = function(x) 2 * drop(w1 * (x - y)))


###################################################
### chunk number 14: 
###################################################
set.seed(12345)
wvec <- 1:9
wmat <- crossprod(matrix(rnorm(81), 9, 9))/9
fit.wls <- activeSet(Atot, "LS", y = y, weights = wvec)
fit.gls <- activeSet(Atot, "GLS", y = y, weights = wmat)


###################################################
### chunk number 15: 
###################################################
fit.qua <- activeSet(Atot, "quantile", y = y, weights = wvec, aw = 0.3, bw = 0.7)


###################################################
### chunk number 16: 
###################################################
fit.abs <- activeSet(Atot, "L1", y = y, weights = w1)


###################################################
### chunk number 17: 
###################################################
fit.eps <- activeSet(Atot, "L1eps", y = y, weights = w1, eps = 1e-04)
fit.pow <- activeSet(Atot, "Lp", y = y, weights = w1, p = 1.2)


###################################################
### chunk number 18: 
###################################################
fit.che <- activeSet(Atot, "chebyshev", y = y, weights = w1)


###################################################
### chunk number 19: 
###################################################
fit.asy <- activeSet(Atot, "asyLS", y = y, weights = w1, aw = 2, bw = 1)


###################################################
### chunk number 20: 
###################################################
fit.hub <- activeSet(Atot, "huber", y = y, weights = w1, eps = 1)
fit.svm <- activeSet(Atot, "SILF", y = y, weights = w1, beta = 0.8, eps = 0.2)


###################################################
### chunk number 21: 
###################################################
set.seed(12345)
yp <- rpois(9, 5)
x0 <- 1:9
fit.poi <- activeSet(Atot, "poisson", x0 = x0, y = yp)


###################################################
### chunk number 22: 
###################################################
Atree <- matrix(c(1, 1, 2, 2, 2, 3, 3, 8, 2, 3, 4, 5, 6, 7, 8, 9), 8, 2)
Atree
fit.tree <- activeSet(Atree, "LS", y = y, weights = w1)


###################################################
### chunk number 23: 
###################################################
Aloop <- matrix(c(1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 3, 3, 4, 5, 6, 6,
  7, 8, 9, 9), 10, 2)
Aloop
fit.loop <- activeSet(Aloop, "LS", y = y, weights = w1)


###################################################
### chunk number 24: 
###################################################
Ablock <- cbind(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3),
  rep(6, 3)), c(rep(c(4, 5, 6), 3), rep(c(7, 8, 9), 3)))
Ablock
fit.block <- activeSet(Ablock, "LS", y = y, weights = w1)


###################################################
### chunk number 25: 
###################################################
pava.fitted <- gpava(1:9, y)$x
aset.fitted <- activeSet(Atot, "LS", weights = w1, y = y)$x
mse <- mean((pava.fitted - aset.fitted)^2)
mse


