library(cvxr)
library(isotone)
data("pituitary", package = "isotone")
head(pituitary)

n <- nrow(pituitary)
y <- matrix(pituitary$size, ncol = 1)
x <- Variable(n)
objective <- Minimize(Pnorm(y - x, 2))
## M <- toeplitz(c(1, -1, rep(0, n-2)))
## M[lower.tri(M)] <- 0
## M[n, n] <- 0
## constraint <- list(M %*% x <= 0)
constraint <- list(diff(x) >= 0)

problem <- Problem(objective, constraint)
result1 <- solve(problem)
sol_x <- result1$getValue(x)

###################################################
### chunk number 3:
###################################################
res1 <- with(pituitary, gpava(age, size, ties = "primary"))
res2 <- with(pituitary, gpava(age, size, ties = "secondary"))
res3 <- with(pituitary, gpava(age, size, ties = "tertiary"))

cbind(res1$x, res2$x, res3$x, sol_x)

## For the primary method we can have different (monotonic) fitted
## values within tied observa- tions, for the secondary the fitted
## values are the same within ties, and the tertiary approach only
## requires monotonicity on the means.

## Secondary

x2 <- Variable(n)
objective2 <- Minimize(Pnorm(y - x2, 2))
secondary_constraints <- lapply(base::split(x = seq_len(nrow(pituitary)),
                                            f = pituitary$age),
                                function(i) diff(x2[i]) == 0)
constraint2 <- c(diff(x2) >= 0,
                 secondary_constraints)

problem2 <- Problem(objective2, constraint2)
result2 <- solve(problem2)
sol_x2 <- result2$getValue(x2)
cbind(res1$x, res2$x, res3$x, sol_x, sol_x2)

## Tertiary

x3 <- Variable(n)
objective3 <- Minimize(Pnorm(y - x3, 2))
blocks <- base::split(x = seq_len(nrow(pituitary)),
                      f = pituitary$age)
block_means <- lapply(blocks,
                      function(i) {
                          n_i <- length(i)
                          v <- numeric(n)
                          v[i] <- 1.0 / n_i
                          M <- matrix(v, nrow = 1)
                          print(M)
                          M %*% x3
                      })
block_mean_vector <- do.call(VStack, block_means)
constraint3 <- list(diff(block_mean_vector) >= 0)

problem3 <- Problem(objective3, constraint3)
result3 <- solve(problem3)
sol_x3 <- result3$getValue(x3)
cbind(res1$x, res2$x, res3$x, sol_x, sol_x2, sol_x3)


gp1 <- gpava(pituitary[,1],pituitary[,2], solver = weighted.mean)
gp2 <- gpava(pituitary[,1],pituitary[,2], solver = weighted.median)
gp3 <- gpava(pituitary[,1],pituitary[,2], solver = weighted.fractile, p = 0.25)

x4 <- Variable(n)
objective4 <- Minimize(Pnorm(y - x4, 1))
constraint4 <- list(diff(x4) >= 0)

problem4 <- Problem(objective4, constraint4)
result4 <- solve(problem4)
sol_x4 <- result4$getValue(x4)
cbind(gp2$x, sol_x4)

### Why it fails?
x5 <- Variable(n)
objective5 <- Minimize(Pnorm(y - x5, 0.25))
constraint5 <- list(diff(x5) >= 0)
problem5 <- Problem(objective5, constraint5)
result5 <- solve(problem5)
sol_x5<- result5$getValue(x5)
cbind(gp3$x, sol_x5)


### Global warming data
library(readr)
d <- read_table(file = "global.txt", skip=16)[, c("YEAR", "ANNUAL")]
n <- nrow(d)
y <- d$ANNUAL
x <- Variable(n)
objective <- Minimize(norm(y - x, 2))
constraint <- list(diff(x) >= 0)
prob <- Problem(objective, constraint)
res <- solve(prob)
fit <- res$getValue(x)

get_iso_estimate <- function(d, indices) {
    boot_data <- d[indices, ]
    boot_data <- boot_data[order(boot_data$YEAR), ]
    y <- boot_data$ANNUAL
    x <- Variable(n)
    objective <- Minimize(norm(y - x, 2))
    constraint <- list(diff(x) >= 0)
    prob <- Problem(objective, constraint)
    res.temp <- solve(prob)
    res.temp$getValue(x)
}
library(boot)
ans <- boot(data = d, statistic = get_iso_estimate, R = 999)
boot.se <- apply(ans$t, 2, sd)
plot.data <- data.frame(year = d$YEAR, annual = d$ANNUAL, est = fit,
                        lower = fit - 1.96 * boot.se, upper = fit + 1.96 * boot.se)
library(ggplot2)
(plot1 <- ggplot(data = plot.data) +
     geom_point(mapping = aes(year, annual), color = "red") +
     geom_line(mapping = aes(year, est), color = "blue") +
     geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3))

n <- nrow(d)
y <- d$ANNUAL
x <- Variable(n)
objective <- Minimize(norm(y - x, 2))
constraint <- list(diff(diff(x)) >= 0)
prob <- Problem(objective, constraint)
res2 <- solve(prob)
fit2 <- res2$getValue(x)

get_iso_estimate2 <- function(d, indices) {
    boot_data <- d[indices, ]
    boot_data <- boot_data[order(boot_data$YEAR), ]
    y <- boot_data$ANNUAL
    x <- Variable(n)
    objective <- Minimize(norm(y - x, 2))
    constraint <- list(diff(x) >= 0)
    prob <- Problem(objective, constraint)
    res.temp <- solve(prob)
    res.temp$getValue(x)
}
ans <- boot(data = d, statistic = get_iso_estimate2, R = 999)
boot.se2 <- apply(ans$t, 2, sd)
plot.data2 <- data.frame(year = d$YEAR, annual = d$ANNUAL, est = fit2,
                         lower = fit2 - 1.96 * boot.se2, upper = fit2 + 1.96 * boot.se2)

(plot2 <- ggplot(data = plot.data2) +
     geom_point(mapping = aes(year, annual), color = "red") +
     geom_line(mapping = aes(year, est), color = "blue") +
     geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3))

##
## Now fit Ryan's model
##

n <- nrow(d)
y <- d$ANNUAL
x <- Variable(n)

## lam <- 0.44
## objective <- Minimize(norm(y - x, 2) + lam * sum(pos(diff(diff(x)))))
## prob3 <- Problem(objective)
## res3 <- solve(prob3)
## fit3 <- res3$getValue(x)

## Nearly isotonic example from the paper

library(readr)
d <- read_table(file = "global.txt", skip=16)[, c("YEAR", "ANNUAL")]
n <- nrow(d)
y <- d$ANNUAL
x <- Variable(n)

get_neariso_estimate <- function(d, indices, lambda) {
    boot_data <- d[indices, ]
    boot_data <- boot_data[order(boot_data$YEAR), ]
    y <- boot_data$ANNUAL
    x <- Variable(n)
    objective <- Minimize(0.5 * norm(y - x, 2) + lambda * sum(pos(diff(x))))
    prob <- Problem(objective)
    solve(prob)$getValue(x)
}

library(boot)
set.seed(129)
## Ryan uses lambda 0.44
system.time(ans <- boot(data = d, statistic = get_neariso_estimate, R = 999, lambda = 0.44))

ci.neariso <- t(sapply(seq_len(n),
                       function(i) boot.ci(boot.out = ans, conf = 0.95,
                                           type = "norm", index = i)$normal[-1]))

data.neariso <- data.frame(year = d$YEAR, annual = d$ANNUAL, est = ans$t0,
                           lower = ci.neariso[, 1], upper = ci.neariso[, 2])

(plot.neariso <- ggplot(data = data.neariso) +
     geom_point(mapping = aes(year, annual), color = "red") +
     geom_line(mapping = aes(year, est), color = "blue") +
     geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3) +
     labs(x = "Year", y = "Temperature Anomalies")
)
pdf("neariso.pdf")
plot.neariso
dev.off()


### REPLICATION
library(isotone)
data("pituitary")
res1 <- with(pituitary, gpava(age, size, ties = "primary"))
res2 <- with(pituitary, gpava(age, size, ties = "secondary"))
res3 <- with(pituitary, gpava(age, size, ties = "tertiary"))

library(cvxr)
## Primary
x_p <- with(pituitary, {
    x <- Variable(length(size))
    objective <- Minimize(Pnorm(size - x, 2))
    constraint <- list(diff(x) >= 0)
    problem <- Problem(objective, constraint)
    result <- solve(problem)
    result$getValue(x)
})

## Secondary
x_s <- with(pituitary, {
    n <- length(size); x <- Variable(n)
    objective <- Minimize(Pnorm(size - x, 2))
    secondary_constraints <- lapply(base::split(x = seq_len(n),
                                                f = age),
                                    function(i) diff(x[i]) == 0)
    constraint <- c(diff(x) >= 0,
                    secondary_constraints)
    problem <- Problem(objective, constraint)
    solve(problem)$getValue(x)
})

## Tertiary
x_t <- with(pituitary, {
    n <- length(size); x <- Variable(n)
    objective <- Minimize(Pnorm(size - x, 2))
    blocks <- base::split(x = seq_len(n),
                          f = pituitary$age)
    block_means <- lapply(blocks, function(i) {
        v <- numeric(n)
        v[i] <- 1.0 / length(i)
        matrix(v, nrow = 1) %*% x
    })
    block_mean_vector <- do.call(VStack, block_means)
    constraint <- list(diff(block_mean_vector) >= 0)
    problem <- Problem(objective, constraint)
    solve(problem)$getValue(x)
})
