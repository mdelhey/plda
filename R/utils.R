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


linear.discriminant <- function(x, mu, sigma, pi) {    
    t(x) %*% solve(sigma) %*% mu - 0.5*t(mu) %*% solve(sigma) %*% mu + log(pi)
}

quadratic.discriminant <- function(x, mu, sigma, pi) {
    (-1/2)*log(det(sigma)) - (-1/2)*t(x - mu) %*% solve(sigma) %*% (x - mu) + log(pi)
}
    
normal.density <- function(x, mu, sigma) {
    det(2*pi*sigma)^(-1/2) * exp(-1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu))
}

#' @title Estimate densities. 
#' @inheritParams plda
estimate.densities <- function(X, mu.hat, sigma.hat, method = "linear") {
    K <- nrow(mu.hat)
    sigma.hat.use <- switch(
        method
      , linear = sigma.hat
      , quadratic = sigma.hat[[k]]
    )
    
    densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
        normal.density(as.matrix(x), mu.hat[k, ], sigma.hat.use)))))
    
    if (any(densities >= 1))
        warning("some estimated densities >= 1!")
}
