library(cvxr)
set.seed(43)
n=100
k=4
lam=.1
alpha=.5
x=rnorm(n)
z=matrix(rnorm(n*3),n,3)
y=2*x+2*x*z[,2]+ .1*rnorm(n)

beta=Variable(4)
xz=matrix(x,nrow(z),ncol(z))*z
objective <- Minimize( (.5/n)*sum((y - x* beta[1] -(xz)%*%beta[2:4] )^2)
+ (1-alpha)*lam*(norm(beta[1:4])+norm(beta[2:4]))+alpha*lam*sum(abs(beta[2:4])))

prob <- Problem(objective)

system.time(result <- solve(prob))

b=result$getValue(beta)


beta1 <- Variable(1)
betaRest <- Variable(3)
betafull <- VStack(beta1, betaRest)
xz=matrix(x,nrow(z),ncol(z))*z
objective <- Minimize( (.5/n)*sum((y - x* beta1 -(xz)%*%betaRest )^2)
+ (1-alpha)*lam*(norm(betafull) + norm(betaRest))+alpha*lam*sum(abs(betaRest)))

prob2 <- Problem(objective)

system.time(result <- solve(prob2))

bb=result$getValue(betafull)
