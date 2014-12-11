devtools::install("~/plda")
devtools::load_all("~/plda")
library(plda)

# Load data
data(iris)

# Data matricies
X <- as.matrix(iris[, c(1,2,3,4)])
y <- as.factor(as.vector(iris[, 5]))

# Train, test data matricies
ind <- split.data(X, train = 0.7, val = 0, test = 0.3)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

# plda::plda
fit.plda <- plda(X.trn, y.trn, method = "linear")
fit.pqda <- plda(X.trn, y.trn, method = "quadratic")
pred.plda <- predict(fit.plda, y.trn)
pred.pqda <- fitted(fit.pqda)

# MASS::lda
fit.lda.formula <- MASS::lda(Species ~ ., data = iris, method = "mle") # Should be the same as bellow
fit.lda <- MASS::lda(X.trn, y, method = "mle")
if (any(coef(fit.lda.formula) != coef(fit.lda)))
    stop("Formula MASS::lda disagrees with matrix MASS::lda")
pred.lda <- predict(fit.lda, X.trn)$class

# e1071::naivebayes
fit.nb <- e1071::naiveBayes(X.trn, y, lapace = 0)
pred.nb <- predict(fit.nb, X.trn)

# PoiClaClu::Classify
fit.poicla <- PoiClaClu::Classify(X.trn, y, xte = X.trn, rho = 0, type = "mle")
pred.poicla <- fit.poicla$ytehat

# Compare results.
which(pred.plda != pred.pqda)
which(pred.plda != pred.lda)
which(pred.plda != pred.nb)
which(pred.plda != pred.poicla)


