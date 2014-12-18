library(plda)
set.seed(1)
################################################################################
### DATA
################################################################################
# Load data
data(iris)

# Data matricies
X <- as.matrix(iris[, c(1,2,3,4)])
y <- as.factor(as.vector(iris[, 5]))

iris.results <- t(replicate(100, {
    # Train, test data matricies
    ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
    X.trn <- X[ind$trn, ]
    y.trn <- y[ind$trn]
    X.tst <- X[ind$tst, ]
    y.tst <- y[ind$tst]

    ################################################################################
    ### MODELS
    ################################################################################
    fit.plda <- plda(X.trn, y.trn, type = "linear",      prior = "uniform")
    fit.pqda <- plda(X.trn, y.trn, type = "quadratic",   prior = "uniform")
    fit.ppda <- plda(X.trn, y.trn, type = "poisson",     prior = "uniform")
    fit.spda <- plda(X.trn, y.trn, type = "poisson.seq", prior = "uniform", size.factor = "quantile")

    pred.plda.trn <- fitted(fit.plda)
    pred.pqda.trn <- fitted(fit.pqda)
    pred.ppda.trn <- fitted(fit.ppda)
    pred.spda.trn <- fitted(fit.spda)

    pred.plda.tst <- predict(fit.plda, X.tst)
    pred.pqda.tst <- predict(fit.pqda, X.tst)
    pred.ppda.tst <- predict(fit.ppda, X.tst)
    pred.spda.tst <- predict(fit.spda, X.tst)

    ################################################################################
    ### MISCLASSIFICATION RATE
    ################################################################################
    errors <- c(
        (1 - length(which(pred.plda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.pqda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.ppda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(pred.spda.tst == y.tst)) / length(y.tst)) * 100
    )
}))

xtable::xtable(cbind(
    apply(iris.results, 2, mean)
  , apply(iris.results, 2, sd)
))

