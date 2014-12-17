# For development
devtools::install("~/plda")
devtools::load_all("~/plda")

set.seed(1)
################################################################################
### Small simulation (for visualization)
################################################################################
n.small <- 50
p.small <- 2
K.small <- 3
lambda.small <- c(10, 28, 48)
sim.small <- simulate.poisson(n = n.small, p = p.small, K = K.small, lambda = lambda.small)

# Data matricies
X <- as.matrix(sim.small$X)
y <- as.factor(sim.small$y)

fit <- plda(X, y, type = "linear", prior = "uniform")

plot(fit, X, y)

# Errors
length(which(fitted(fit) == y)) / length(y)

ggplot() +
  geom_point(aes(x = X[, 1], y = X[, 2], color = fitted(fit))) +
  theme_bw() + xlab("feature 1") + ylab("feature 2") + 
  xlim(Xmin, Xmax) + ylim(Ymin, Ymax) 


################################################################################
### Large simulation (for accuracy)
################################################################################

n.large <- 100
p.large <- 25
K.large <- 5
lambda.large <- c(10, 12, 28, 48, 100)
sim.small <- simulate.poisson(n = n.large, p = p.large, K = K.large, lambda = lambda.large)

# Data matricies
X <- as.matrix(sim.small$X)
y <- as.factor(sim.small$y)

# Train, test data matricies
ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

# Fit models
fit.lda <- plda(X.trn, y.trn, type = "linear",    prior = "uniform") 
fit.qda <- plda(X.trn, y.trn, type = "quadratic", prior = "uniform")
fit.pda <- plda(X.trn, y.trn, type = "poisson",   prior = "uniform", size.factor = "mle")
fit.poicla <- PoiClaClu::Classify(x = X.trn, y = y.trn, xte = X.tst)

# Classify test data
pred.lda <- predict(fit.lda, X.tst)
pred.qda <- predict(fit.qda, X.tst)
pred.pda <- predict(fit.pda, X.tst)
pred.poicla <- fit.poicla$ytehat

# Errors
length(which(pred.lda == y.tst)) / length(y.tst)
length(which(pred.qda == y.tst)) / length(y.tst)
length(which(pred.pda == y.tst)) / length(y.tst)
length(which(pred.poicla == y.tst)) / length(y.tst)

################################################################################
### Compare with PoiClaClu
################################################################################
# Compare with PoiClaClu
#cv <- PoiClaClu::Classify.cv(X, y)
fit.poicla <- PoiClaClu::Classify(X, y, xte = X.grid, type = "mle", transform = FALSE, rho = 0)
#pred.trn <- PoiClaClu::Classify(X, y, xte = X, type = "mle", transform = TRUE, rho = 0)$ytehat
y.grid <- fit.poicla$ytehat

sim <- PoiClaClu::CountDataSet(n = 50, p = 1000, sdsignal = 0.025, K = 3, param = 0.1)

X <- as.matrix(sim$x)
y <- as.factor(sim$y)

ind <- split.data(X, train = 0.5, val = 0, test = 0.5)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

# Fit models
fit.lda <- plda(X.trn, y.trn, type = "linear",    prior = "uniform") 
#fit.qda <- plda(X.trn, y.trn, type = "quadratic", prior = "uniform")
#fit.pda <- plda(X.trn, y.trn, type = "poisson",   prior = "uniform", size.factor = "mle")
fit.poicla <- PoiClaClu::Classify(X.trn, y.trn, xte = X.tst)

# Classify test data
pred.lda <- predict(fit.lda, X.tst)
#pred.qda <- predict(fit.qda, X.tst)
#pred.pda <- predict(fit.pda, X.tst)
pred.poicla <- fit.poicla$ytehat

length(which((pred.lda) != y.tst)) / length(y.tst)
length(which((pred.poicla) != y.tst)) / length(y.tst)
