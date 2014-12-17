# For development
devtools::install("~/plda")
devtools::load_all("~/plda")

set.seed(1)
################################################################################
### Splice Junction Gene sequences
################################################################################
fname <- system.file("extdata/splice-junction-gene-sequences/splice.data", package = "plda")
splice <- read.csv(fname, header = FALSE, stringsAsFactors = FALSE)
names(splice) <- c("class", "instance", "sequence")

X <- count.sequence(splice$sequence, words = 1)
y <- as.factor(splice$class)

# Train, test data matricies
ind <- split.data(X, train = 0.5, val = 0, test = 0.5)
X.trn <- X[ind$trn, ]
y.trn <- y[ind$trn]
X.tst <- X[ind$tst, ]
y.tst <- y[ind$tst]

fit.lda <- plda(X.trn, y.trn, type = "linear")
fit.qda <- plda(X.trn, y.trn, type = "quadratic")
fit.pda <- plda(X.trn, y.trn, type = "poisson")

y.lda.tst <- predict(fit.lda, X.tst)
y.qda.tst <- predict(fit.qda, X.tst)
y.pda.tst <- predict(fit.pda, X.tst)

# Accuracy
length(which(y.lda.tst == y.tst)) / length(y.tst)
length(which(y.qda.tst == y.tst)) / length(y.tst)
length(which(y.pda.tst == y.tst)) / length(y.tst)
