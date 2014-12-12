#' @title Poisson discriminant analysis.
#' @param X Matrix of predictors.
#' @param y Factor vector of class or group labels.
#' @param type Type of discriminant analysis to be performed. Linear assumes that all population
#' class covariances are equal. Quadratic estimates each population class covariance
#' separately. Linear and Quadratic assume predictors are normally distributed (i.e. continous
#' data). Poisson models count data.
#' @return List of class "plda" containing the classifier results.
plda <- function(X, y, type = c("linear", "quadratic", "poisson"), prior = c("proportion", "uniform")) {
    stopifnot(is.matrix(X), is.factor(y), length(y) == nrow(X))
    type  <- match.arg(type)
    prior <- match.arg(prior)
    
    # House-keeping variables.
    K <- length(levels(y))
    N.k <- t(as.matrix(by(X, y, nrow)))
    N <- nrow(X)
    P <- ncol(X)
    beta <- 0 # Smoothing parameter for Poisson. 0 is MLE

    # Estimate prior. 
    pi.hat <- switch(
        prior
      , uniform    = setNames(rep(1 / K, K), levels(y))
      , proportion = N.k / N
    )
    
    if (type == "poisson") {
        # Dot notation
        Xi. <- apply(X, 1, sum)
        X.j <- apply(X, 2, sum)
        X.. <- sum(X)

        # MLE Estimates
        s.hat <- Xi. / X.. # size factor
        g.hat <- X.j
        N.hat <- as.matrix(s.hat) %*% t(as.matrix(g.hat))
        a <- as.matrix(do.call("rbind", by(X, y, function(x) colSums(x) + beta)))
        b <- as.matrix(do.call("rbind", by(X, y, function(x) colSums(N.hat) + beta)))
        d.hat <- a / b # dkj > 1, jth feature is over-expressed relative to kth class       
        
        if (!all.equal(sum(s.hat), 1))
            stop("Error estimating s.hat, doesn't sum to 1.")        

        # Calculate Liklihood of Poisson using row-class dependent lambda
        densities <- estimate.poisson.densities(X, N.hat, d.hat)

        # Use bayes rule to calculate posterior densities
        posteriors <- bayes.rule(X, densities, pi.hat)
        fitted.posteriors <- factor(levels(y)[apply(posteriors, 1, which.max)], levels = levels(y))

        # Return poisson estimates
        estimates = list(
            pi.hat = pi.hat
          , N.hat  = N.hat
          , d.hat  = d.hat
        )        
    }
    
    if (type %in% c("linear", "quadratic")) {
        # MLE Estimates
        mu.hat <- as.matrix(do.call("rbind", by(X, y, colMeans)))
        sigma.hat <- switch(
            type
          , linear = Reduce("+", by(X, y, cov)) / K
          , quadratic = by(X, y, cov)
        )

        # Calculate Liklihood of MVN
        densities <- estimate.normal.densities(X, mu.hat, sigma.hat, type = type)

        # Use bayes rule to calculate posterior densities
        posteriors <- bayes.rule(X, densities, pi.hat)        
        fitted.posteriors <- factor(levels(y)[apply(posteriors, 1, which.max)], levels = levels(y))

        # Compare with Fisher's discriminant rule to assign classes
        if (type == "linear") 
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                linear.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat, pi.hat[k])))))
        if (type == "quadratic")
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                quadratic.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat[[k]], pi.hat[k])))))
        
        # Ensure we make the same predictions.
        fitted.discriminants <- factor(levels(y)[apply(discriminants, 1, which.max)], levels = levels(y))        
        
        if (!isTRUE(all.equal(fitted.discriminants, fitted.posteriors)))
            warning("discrimants and posteriors do not agree!")

        # Return normal-dist estimates
        estimates = list(
            pi.hat    = pi.hat
          , mu.hat    = mu.hat
          , sigma.hat = sigma.hat
        )
    }
    
    result <- list(
        parameters = list(
            P = P
          , N = N
          , K = K
          , N.k = N.k
          , type = type
          , levels = levels(y)
        )
      , estimates = estimates
      , fitted = fitted.posteriors
      , posteriors = posteriors
      , densities = densities        
    )
    class(result) <- "plda"
    return(result)
}
