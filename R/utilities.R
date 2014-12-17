#' @title Split data.
#' @param x Matrix or vector to split.
#' @param train Proportion of rows to label as train.
#' @param val Proportion of rows to label as validation.
#' @param test Proportion of rows to label as test.
#' @return List of split indicies.
split.data <- function(x, train = 0.7, val = 0.0, test = 0.3) {
    if (!(is.matrix(x) || is.vector(x) || is.data.frame(x)))
        stop("X must be a matrix, dataframe, or vector.")
    if (1 != train + val + test)
        stop("Train, validation, and test proportions must sum to 1.")

    n <- nrow(x)
    n.trn <- floor(train * n)
    n.val <- floor(val * n)
    n.tst <- n - (n.trn + n.val)
    
    ind.trn <- sample(n, size = n.trn, replace = FALSE)
    ind.val <- sample(n, size = n.val, replace = FALSE)
    ind.tst <- sample(n, size = n.tst, replace = FALSE)

    return(list(trn = ind.trn, val = ind.val, tst = ind.tst))
}

#' @title Simulate Poisson.
#' @return List of X and y.
#' X is (n*K by p) data matrix composed of independent univariate poissons.
#' y is (n by 1) factor of class labels.
simulate.poisson <- function(n, p, K, lambda) {
    stopifnot(is.numeric(n), is.vector(n))
    stopifnot(is.numeric(p), is.vector(p))
    stopifnot(is.numeric(K), is.vector(K))
    if (length(lambda) != K || !is.vector(lambda))
        stop("lambda must be a vector of length K")    

    X <- do.call("rbind", lapply(1:K, function(k)
        replicate(p, rpois(n, lambda = lambda[k]))))
    
    y <- as.factor(do.call("c", lapply(1:K, function(k)
        rep(k, n))))
    
    simulation <- list(X = X, y = y, K = K, lambda = lambda)
    
    return(simulation)
}

#' @title Point grid.
point.grid <- function(Xmin, Xmax, Ymin, Ymax, n.seq = 1000) {
    x <- seq(Xmin, Xmax, length.out = n.seq)
    y <- seq(Ymin, Ymax, length.out = n.seq)
    as.matrix(expand.grid(x, y))
}

#' @title Geometric mean.
geo.mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

count.sequence <- function(sequence, words = 2) {
    seqs <- strsplit(tolower(stringr::str_trim(sequence)), "")
    counts <- do.call("rbind", (lapply(1:length(seqs), function(i)
        seqinr::count(seqs[[i]], words))))

    return(counts)
}
