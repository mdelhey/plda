path <- path.package("plda")
devtools::install(path)
devtools::load_all(path)
library(plda)

library(ggplot2)
set.seed(1)

################################################################################
### Small simulation (for visualization)
################################################################################
#sim <- PoiClaClu::CountDataSet(n = 100, p = 2, sdsignal = 0.5, K = 3, param = 10)
n.small <- 50
p.small <- 2
K.small <- 2
lambda.small <- c(10, 28)
# X = (100 x 2), K = 2
sim.small <- simulate.poisson(n = n.small, p = p.small, K = K.small, lambda = lambda.small)

# Data matricies
X <- as.matrix(sim$X)
y <- as.factor(sim$y)

ggplot() + geom_point(aes(x = X[, 1], y = X[, 2], color = y)) +
    theme_bw() + xlab("feature 1") + ylab("feature 2") + 
        xlim(0, 40) + ylim(0, 45)

Xmin <- 0
Xmax <- 40



################################################################################
### Large simulation (for accuracy)
################################################################################

# Train, test data matricies
ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

