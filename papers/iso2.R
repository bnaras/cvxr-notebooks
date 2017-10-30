## Nearly convex fit from paper
library(cvxr)
library(readr)
d <- read_table(file = "global.txt", skip=16)[, c("YEAR", "ANNUAL")]
n <- nrow(d)
y <- d$ANNUAL
x <- Variable(n)


get_nearconvex_estimate <- function(d, indices, lambda) {
    boot_data <- d[indices, ]
    boot_data <- boot_data[order(boot_data$YEAR), ]
    y <- boot_data$ANNUAL
    x <- Variable(n)
    objective <- Minimize(0.5 * norm(y - x, 2) + lambda * sum(pos(diff(diff(x)))))
    prob <- Problem(objective)
    solve(prob)$getValue(x)
}

get_nearconvex_estimate2 <- function(d, indices, lambda) {
    boot_data <- d[indices, ]
    boot_data <- boot_data[order(boot_data$YEAR), ]
    y <- boot_data$ANNUAL
    x <- Variable(n)
    objective <- Minimize(0.5 * norm(y - x, 2) + lambda * sum(pos(diff(x, differences = 2))))
    prob <- Problem(objective)
    solve(prob)$getValue(x)
}

system.time(res1 <- get_nearconvex_estimate(d, seq_len(nrow(d)), lambda=0.44))
##   user  system elapsed
##349.975   1.335 352.009

system.time(res2 <- get_nearconvex_estimate2(d, seq_len(nrow(d)), lambda=0.44))
##   user  system elapsed
##351.696   1.137 353.308

library(boot)
set.seed(2829)
system.time(ans <- boot(data = d, statistic = get_nearconvex_estimate, R = 100, lambda = 0.44))

ci.nearconvex <- t(sapply(seq_len(n),
                          function(i) boot.ci(boot.out = ans, conf = 0.95,
                                              type = "norm", index = i)$normal[-1]))

data.nearconvex <- data.frame(year = d$YEAR, annual = d$ANNUAL, est = ans$t0,
                              lower = ci.nearconvex[, 1], upper = ci.nearconvex[, 2])

library(ggplot2)
(plot.nearconvex <- ggplot(data = data.nearconvex) +
     geom_point(mapping = aes(year, annual), color = "red") +
     geom_line(mapping = aes(year, est), color = "blue") +
     geom_ribbon(mapping = aes(x = year, ymin = lower,ymax = upper),alpha=0.3) +
     labs(x = "Year", y = "Temperature Anomalies")
)

pdf("nearconvex.pdf")
plot.nearconvex
dev.off()

