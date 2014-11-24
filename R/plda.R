plda <- function(X, y, method = c("linear", "quadratic"), type = c("normal", "poisson")) {
    stopifnot(is.matrix(X), is.factor(y), length(y) == nrow(X))
    method <- arg.

    # House-keeping variables.
    K <- length(levels(y))
    N.k <- t(as.matrix(by(X, y, nrow)))
    N <- nrow(X)
    P <- ncol(X)
    
    # MLE Estimates
    pi.hat <- N.k / N
    mu.hat <- as.matrix(do.call("rbind", by(X, y, colMeans)))
    sigma.hat <- switch(
        method
      , linear = Reduce("+", by(X, y, cov)) / K
      , quadratic = by(X, y, cov)
    )    
    
    #diag(1/sqrt(diag(var(X - mu.hat[y, ]))))     # Coefficients?

    # Calculate Liklihood of MVN
    if (method == "linear")
        densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
            normal.density(as.matrix(x), mu.hat[k, ], sigma.hat)))))
    if (method == "quadratic")
        densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
            normal.density(as.matrix(x), mu.hat[k, ], sigma.hat[[k]])))))

    # Use bayes rule to calculate posterior densities
    posteriors <- do.call("cbind", lapply(1:K, function(j)
        unlist(lapply(1:N, function(i)    
            (pi.hat[j] * densities[i, j]) / sum(pi.hat * densities[i, ])))))

    if (!isTRUE(all.equal(apply(posteriors, 1, sum), rep(1, N), check.attributes = FALSE)))
        warning("posteriors do not sum to one!")
    
    fitted.posteriors <- factor(levels(y)[apply(posteriors, 1, which.max)], levels = levels(y))

    # Compare with Fisher's discriminant rule to assign classes
    if (method == "linear") {
        discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
            linear.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat, pi.hat[k])))))
        
        # Ensure we make the same predictions.
        fitted.discriminants <- factor(levels(y)[apply(discriminants, 1, which.max)], levels = levels(y))

        if (!(all.equal(fitted.discriminants, fitted.posteriors)))
            warning("discrimants and posteriors do not agree!")
    }    
    
    result <- list(
        parameters = list(
            P = P
          , N = N
          , K = K
          , N.k = N.k
          , method = method
        )
      , fitted = fitted.posteriors
      , posteriors = posteriors
      , densities = densities
      , pi.hat = pi.hat
      , mu.hat = mu.hat
      , sigma.hat = sigma.hat
    )
    class(result) <- "plda"
    return(result)
}

predict.plda <- function(result, ...) {
    
}

linear.discriminant <- function(x, mu, sigma, pi) {    
    t(x) %*% solve(sigma) %*% mu - 0.5*t(mu) %*% solve(sigma) %*% mu + log(pi)
}
    
normal.density <- function(x, mu, sigma) {
    det(2*pi*sigma)^(-1/2) * exp(-1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu))
}

print.plda <- function(result, ...) {
    cat(sprintf("X(%i x %i) with K = %i\n", result$parameters$N, result$parameters$P, result$parameters$K))
    cat("Method: %s DA\n", method)

    cat("\npi.hat:\n")
    printCoefmat(t(as.matrix(result$pi.hat)))

    cat("\nmu.hat:\n")
    printCoefmat(result$mu.hat)    
}

