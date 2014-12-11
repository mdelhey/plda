#' @title Bayes rule for estimating posteriors.
#' @inheritParams plda
bayes.rule <- function(X, densities, pi) {
    K <- ncol(densities)
    N <- nrow(X)
    
    posteriors <- do.call("cbind", lapply(1:K, function(j)
        unlist(lapply(1:N, function(i)    
            (pi[j] * densities[i, j]) / sum(pi * densities[i, ])))))

    if (!isTRUE(all.equal(apply(posteriors, 1, sum), rep(1, N), check.attributes = FALSE)))
        stop("posteriors do not sum to one!")

    return(posteriors)
}

#' @title Estimate densities. 
#' @inheritParams plda
estimate.densities <- function(X, mu.hat, sigma.hat, method = "linear") {
    K <- nrow(mu.hat)    
    sigma.hat.use <- switch(
        method
      , linear = "sigma.hat"
      , quadratic = "sigma.hat[[k]]"
    )
    
    densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
        normal.density(as.matrix(x), mu.hat[k, ], eval(parse(text = sigma.hat.use)))))))

    return(densities)
}


linear.discriminant <- function(x, mu, sigma, pi) {    
    t(x) %*% solve(sigma) %*% mu - 0.5*t(mu) %*% solve(sigma) %*% mu + log(pi)
}

quadratic.discriminant <- function(x, mu, sigma, pi) {
    (-1/2)*log(det(sigma)) - (1/2)*t(x - mu) %*% solve(sigma) %*% (x - mu) + log(pi)
}

#' @title Multivariate normal density
normal.density <- function(x, mu, sigma) {
    det(2*pi*sigma)^(-1/2) * exp(-1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu))
}

#' @title Poisson univariate density.
poisson.density <- function(x, lambda) {
    (lambda^x / factorial(x)) * exp(1)^(-lambda)
}

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
