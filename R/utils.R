#' @title Bayes rule for estimating posteriors.
#' @inheritParams plda
#' @param type Calculate posteriors using raw products of probabilities or summation of logs of
#' probabilities? "Sum" may offer superior floating point accuracy.
bayes.rule <- function(X, densities, pi, type = c("prod", "sum")) {
    type <- match.arg(type)
    K <- ncol(densities)
    N <- nrow(X)

    if (type == "prod") {
        posteriors <- do.call("cbind", lapply(1:K, function(j)
            unlist(lapply(1:N, function(i)    
                (pi[j] * densities[i, j]) / sum(pi * densities[i, ])))))
    }

    if (type == "sum") {
        log.posteriors <- do.call("cbind", lapply(1:K, function(j)
            unlist(lapply(1:N, function(i)
                log(pi[j]) + log(densities[i, j]) - log(sum(pi * densities[i, ]))))))
        posteriors <- exp(1)^log.posteriors
    }

    if (!isTRUE(all.equal(apply(posteriors, 1, sum), rep(1, N), check.attributes = FALSE)))
        stop("posteriors do not sum to one!")

    return(posteriors)
}

#' @title Estimate normal densities. 
#' @inheritParams plda
estimate.normal.densities <- function(X, mu.hat, sigma.hat, type = c("linear", "quadratic")) {
    type  <- match.arg(type)
    K <- nrow(mu.hat)
    
    sigma.hat.use <- switch(
        type
      , linear = "sigma.hat"
      , quadratic = "sigma.hat[[k]]"
    )
    
    densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
        normal.density(as.matrix(x), mu.hat[k, ], eval(parse(text = sigma.hat.use)))))))

    return(densities)
}

#' @title Estimate poisson densities.
#' @inheritParams plda
estimate.poisson.densities <- function(X, N.hat, d.hat) {
    K <- nrow(d.hat)
    N <- nrow(X)
    
    densities <- do.call("rbind", lapply(1:N, function(i) unlist(lapply(1:K, function(k)
        prod(poisson.density(X[i, ], lambda = N.hat[i, ] * d.hat[k, ]))))))

    return(densities)
}

#' @title Calculate discriminants
#' @description Wrapper for discriminant functions. Not currently in use.
calculate.discriminants <- function(X, parameters, type = c("linear", "quadratic", "poisson")) {
    type <- match.arg(type)
    K <- ncol(X)
    N <- nrow(X)
    
    discriminants <- do.call("rbind", lapply(1:N, function(i) unlist(lapply(1:K, function(k)
        poisson.discriminant(X[i, ], s.hat[i], g.hat, d.hat[k, ], pi.hat[k])))))
    
    return(discriminants)
}

#' @title Linear discriminant
linear.discriminant <- function(x, mu, sigma, pi) {    
    t(x) %*% solve(sigma) %*% mu - 0.5*t(mu) %*% solve(sigma) %*% mu + log(pi)
}

#' @title Quadratic discriminant
quadratic.discriminant <- function(x, mu, sigma, pi) {
    (-1/2)*log(det(sigma)) - (1/2)*t(x - mu) %*% solve(sigma) %*% (x - mu) + log(pi)
}

#' @title Poisson discriminant 
poisson.discriminant <- function(x, s, g, d, pi) {
    sum(x*log(d)) - s*sum(g*d) + log(pi)
}

#' @title Multivariate normal density
normal.density <- function(x, mu, sigma) {
    det(2*pi*sigma)^(-1/2) * exp(-1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu))
}

#' @title Poisson univariate density.
poisson.density <- function(x, lambda) {
    (lambda^x / factorial(x)) * exp(1)^(-lambda)
}

fit.posteriors <- function(posteriors, levels) {
    stopifnot(is.character(levels), is.vector(levels), is.matrix(posteriors))
    factor(levels[apply(posteriors, 1, which.max)], levels = levels)
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

point.grid <- function(Xmin, Xmax, Ymin, Ymax, n.seq = 1000) {
    x <- seq(Xmin, Xmax, length.out = n.seq)
    y <- seq(Ymin, Ymax, length.out = n.seq)

    z <- do.call("cbind", lapply(1:length(x), function(i)
        cbind(x[i], y)))
           
    X <- as.matrix()
    z <- expand.grid(x, y)
}
