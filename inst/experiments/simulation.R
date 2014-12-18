library(plda)
set.seed(1)
################################################################################
### Small simulation (for visualization)
################################################################################
n.small <- 50
p.small <- 2
K.small <- 3
lambda.small <- c(10, 28, 48)
sim.small <- simulate.poisson(n = n.small, p = p.small, K = K.small, lambda = lambda.small)

X <- as.matrix(sim.small$X)
y <- as.factor(sim.small$y)

fit.lda <- plda(X, y, type = "linear",      prior = "uniform")
fit.qda <- plda(X, y, type = "quadratic",   prior = "uniform")
fit.pda <- plda(X, y, type = "poisson",     prior = "uniform")
fit.sda <- plda(X, y, type = "poisson.seq", prior = "uniform", size.factor = "medratio")

errors <- c(
    1 - length(which(fitted(fit.lda) == y)) / length(y)
  , 1 - length(which(fitted(fit.qda) == y)) / length(y)
  , 1 - length(which(fitted(fit.pda) == y)) / length(y)
  , 1 - length(which(fitted(fit.sda) == y)) / length(y)    
)    

# Plots
plot.size <- c(4, 4) # width x height

plot(fit.lda, X, y) + ggplot2::ggtitle("linear")
ggplot2::ggsave("decision-boundary-lda.pdf", width = plot.size[1], height = plot.size[2])

plot(fit.qda, X, y) + ggplot2::ggtitle("quadratic")
ggplot2::ggsave("decision-boundary-qda.pdf", width = plot.size[1], height = plot.size[2])

plot(fit.pda, X, y) + ggplot2::ggtitle("poisson")
ggplot2::ggsave("decision-boundary-pda.pdf", width = plot.size[1], height = plot.size[2])

plot(fit.sda, X, y) + ggplot2::ggtitle("poisson sequence")
ggplot2::ggsave("decision-boundary-sda.pdf", width = plot.size[1], height = plot.size[2])

################################################################################
### Large simulation (for accuracy)
################################################################################
sim.large.results <- t(replicate(100, {    
    n.large <- 100
    p.large <- 15
    K.large <- 5
    lambda.large <- c(10, 12, 28, 48, 100)
    sim.large <- simulate.poisson(n = n.large, p = p.large, K = K.large, lambda = lambda.large)

    X <- as.matrix(sim.large$X)
    y <- as.factor(sim.large$y)

    ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
    X.trn <- X[ind$trn, ]
    y.trn <- y[ind$trn]
    X.tst <- X[ind$tst, ]
    y.tst <- y[ind$tst]

    fit.lda <- plda(X.trn, y.trn, type = "linear",      prior = "uniform") 
    fit.qda <- plda(X.trn, y.trn, type = "quadratic",   prior = "uniform")
    fit.pda <- plda(X.trn, y.trn, type = "poisson",     prior = "uniform")
    fit.sda <- plda(X.trn, y.trn, type = "poisson.seq", prior = "uniform", size.factor = "medratio")

    pred.lda <- predict(fit.lda, X.tst)
    pred.qda <- predict(fit.qda, X.tst)
    pred.pda <- predict(fit.pda, X.tst)
    pred.sda <- predict(fit.sda, X.tst)

    errors <- c(
        (1 - length(which(pred.lda == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.qda == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.pda == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.sda == y.tst)) / length(y.tst)) * 100
    )
}))

xtable::xtable(cbind(
    apply(sim.large.results, 2, mean)
  , apply(sim.large.results, 2, sd)
))

