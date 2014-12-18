library(plda)
set.seed(1)
################################################################################
### Splice Junction Gene sequences
################################################################################
fname <- system.file("extdata/splice-junction-gene-sequences/splice.data", package = "plda")
splice <- read.csv(fname, header = FALSE, stringsAsFactors = FALSE)
names(splice) <- c("class", "instance", "sequence")

X <- count.sequence(splice$sequence, words = 2)
y <- as.factor(splice$class)

# Train, test data matricies
seq.results <- t(replicate(10, {
    ind <- split.data(X, train = 0.6, val = 0, test = 0.4)
    X.trn <- X[ind$trn, ]
    y.trn <- y[ind$trn]
    X.tst <- X[ind$tst, ]
    y.tst <- y[ind$tst]

    fit.lda <- plda(X.trn, y.trn, type = "linear",      prior = "proportion")
    fit.qda <- plda(X.trn, y.trn, type = "quadratic",   prior = "proportion")
    fit.pda <- plda(X.trn, y.trn, type = "poisson",     prior = "proportion")
    fit.sda <- plda(X.trn, y.trn, type = "poisson.seq", prior = "proportion", size.factor = "medratio")

    y.lda.tst <- predict(fit.lda, X.tst)
    y.qda.tst <- predict(fit.qda, X.tst)
    y.pda.tst <- predict(fit.pda, X.tst)
    y.sda.tst <- predict(fit.sda, X.tst)

    # Misclassification
    errors <- c(
        (1 - length(which(y.lda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(y.qda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(y.pda.tst == y.tst)) / length(y.tst)) * 100
      , (1 - length(which(y.sda.tst == y.tst)) / length(y.tst)) * 100        
    )
}))

xtable::xtable(cbind(
    apply(seq.results, 2, mean)
  , apply(seq.results, 2, sd)
))
