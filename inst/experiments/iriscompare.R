path <- "~/plda"
devtools::install(path)
devtools::load_all(path)

library(plda)
################################################################################
### DATA
################################################################################
# Load data
data(iris)

# Data matricies
X <- as.matrix(iris[, c(1,2,3,4)])
y <- as.factor(as.vector(iris[, 5]))

# Train, test data matricies
ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

################################################################################
### plda::plda
################################################################################
fit.plda <- plda(X.trn, y.trn, type = "linear",      prior = "uniform")
fit.pqda <- plda(X.trn, y.trn, type = "quadratic",   prior = "uniform")
fit.ppda <- plda(X.trn, y.trn, type = "poisson",     prior = "uniform")
fit.spda <- plda(X.trn, y.trn, type = "poisson.seq", prior = "uniform")

pred.plda.trn <- fitted(fit.plda)
pred.pqda.trn <- fitted(fit.pqda)
pred.ppda.trn <- fitted(fit.ppda)
pred.spda.trn <- fitted(fit.spda)

pred.plda.tst <- predict(fit.plda, X.tst)
pred.pqda.tst <- predict(fit.pqda, X.tst)
pred.ppda.tst <- predict(fit.ppda, X.tst)


################################################################################
### MASS::lda
################################################################################
fit.lda.formula <- MASS::lda(Species ~ ., data = iris, method = "mle") # Should be the same as bellow
fit.lda <- MASS::lda(X.trn, y.trn, method = "mle")
if (any(coef(fit.lda.formula) != coef(fit.lda)))
    stop("Formula MASS::lda disagrees with matrix MASS::lda")
pred.lda.trn <- predict(fit.lda, X.trn)$class
pred.lda.tst <- predict(fit.lda, X.tst)$class

################################################################################
### e1071::naivebayes 
################################################################################
fit.nb <- e1071::naiveBayes(X.trn, y.trn, lapace = 0)
pred.nb.trn <- predict(fit.nb, X.trn)
pred.nb.tst <- predict(fit.nb, X.tst)

################################################################################
### PoiClaClu::Classify
################################################################################
fit.poicla <- PoiClaClu::Classify(X.trn, y.trn, xte = X.trn, type = "mle", transform = FALSE, rho = 0)
pred.poicla.trn <- fit.poicla$ytehat
pred.poicla.tst <- PoiClaClu::Classify(X.trn, y.trn, xte = X.tst, rho = 0, type = "mle")$ytehat

################################################################################
### Compare results
################################################################################
# Training classification
which(pred.plda.trn != pred.pqda.trn)
which(pred.plda.trn != pred.ppda.trn)
which(pred.plda.trn != pred.lda.trn)
which(pred.plda.trn != pred.nb.trn)
which(pred.plda.trn != pred.poicla.trn)

which(pred.ppda.trn != pred.poicla.trn)
which(pred.ppda.tst != pred.poicla.tst)

# Test accuracy
length(which(pred.plda.tst   == y.tst)) / length(y.tst)
length(which(pred.pqda.tst   == y.tst)) / length(y.tst)
length(which(pred.lda.tst    == y.tst)) / length(y.tst)
length(which(pred.poicla.tst == y.tst)) / length(y.tst)

length(which(pred.ppda.tst   == y.tst)) / length(y.tst)
