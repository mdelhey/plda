#' @title Poisson discriminant analysis.
#' @param X Matrix of predictors.
#' @param y Factor vector of class or group labels.
#' @param method Method of discriminant analysis to be
#' performed. Linear assumes that all population class covariances are
#' equal. Quadratic estimates each population class covariance
#' separately.
#' @param dist Underlying distribution of predictors. Normal models
#' continuous data. Poisson models count data.
#' @return List of class "plda" containing the classifier results.
plda <- function(X, y, method = c("linear", "quadratic"), dist = c("normal", "poisson")) {
    stopifnot(is.matrix(X), is.factor(y), length(y) == nrow(X))
    method <- match.arg(method)
    dist <- match.arg(dist)
    
    # House-keeping variables.
    K <- length(levels(y))
    N.k <- t(as.matrix(by(X, y, nrow)))
    N <- nrow(X)
    P <- ncol(X)

    if (dist == "normal") {       
        # MLE Estimates
        pi.hat <- N.k / N
        mu.hat <- as.matrix(do.call("rbind", by(X, y, colMeans)))
        sigma.hat <- switch(
            method
          , linear = Reduce("+", by(X, y, cov)) / K
          , quadratic = by(X, y, cov)
        )

        # Calculate Liklihood of MVN
        densities <- estimate.densities(X, mu.hat, sigma.hat, method = method)
        ## if (method == "linear")
        ##     densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
        ##         normal.density(as.matrix(x), mu.hat[k, ], sigma.hat)))))
        ## if (method == "quadratic")
        ##     densities <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
        ##         normal.density(as.matrix(x), mu.hat[k, ], sigma.hat[[k]])))))

        # Use bayes rule to calculate posterior densities
        posteriors <- bayes.rule(X, densities, pi.hat)            
        fitted.posteriors <- factor(levels(y)[apply(posteriors, 1, which.max)], levels = levels(y))

        # Compare with Fisher's discriminant rule to assign classes
        if (method == "linear") 
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                linear.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat, pi.hat[k])))))
        if (method == "quadratic")
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                quadratic.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat[[k]], pi.hat[k])))))
        
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
          , dist = dist
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
