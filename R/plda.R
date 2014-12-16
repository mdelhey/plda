#' @title Poisson discriminant analysis.
#' @param X Matrix of predictors.
#' @param y Factor vector of class or group labels.
#' @param type Type of discriminant analysis to be performed. Linear assumes that all population
#' class covariances are equal. Quadratic estimates each population class covariance
#' separately. Linear and Quadratic assume predictors are normally distributed (i.e. continous
#' data). Poisson models count data.
#' @return List of class "plda" containing the classifier results.
plda <- function(X, y,
                 type = c("linear", "quadratic", "poisson"),
                 size.factor = c("mle", "quantile", "medratio"),
                 prior = c("uniform", "proportion"))
{
    stopifnot(is.matrix(X), is.factor(y), length(y) == nrow(X))
    type  <- match.arg(type)
    prior <- match.arg(prior)
    size.factor <- match.arg(size.factor)
    
    # House-keeping variables.
    K <- length(levels(y))
    N.k <- t(as.matrix(by(X, y, nrow)))
    N <- nrow(X)
    P <- ncol(X)

    # Estimate prior. 
    pi.hat <- switch(
        prior
      , uniform    = setNames(rep(1 / K, K), levels(y))
      , proportion = N.k / N
    )
    
    if (type == "poisson") {
        # Smoothing parameter for Poisson. 0 is MLE
        beta <- 1

        # Parameter estimates
        s.estimate <- estimate.size.factor(X, type = size.factor)
        s.hat <- s.estimate$s.hat                           # s = counts per obs
        g.hat <- colSums(X)                                 # g = counts per feature                
        N.hat <- matrix(s.hat) %*% matrix(g.hat, nrow = 1)  # N = lambda's
        
        a <- as.matrix(do.call("rbind", by(X, y, function(x) colSums(x) + beta)))
        b <- as.matrix(do.call("rbind", by(N.hat, y, function(x) colSums(x) + beta)))
        d.hat <- a / b # dkj > 1, jth feature is over-expressed relative to kth class       
        
        # Calculate Liklihood of Poisson using row-class dependent lambda
        densities <- estimate.poisson.densities(X, N.hat, d.hat)

        # Use bayes rule to calculate posterior densities
        posteriors <- bayes.rule(X, densities, pi.hat)
        fitted.posteriors <- fit.posteriors(posteriors, levels(y))

        # Compare with discriminant rule to assign classes
        discriminants <- do.call("rbind", lapply(1:N, function(i) unlist(lapply(1:K, function(k)
            poisson.discriminant(X[i, ], s.hat[i], g.hat, d.hat[k, ], pi.hat[k])))))

        # Ensure we make the same predictions.
        fitted.discriminants <- fit.posteriors(discriminants, levels(y))
        
        if (!isTRUE(all.equal(fitted.discriminants, fitted.posteriors)))
            warning("discrimants and posteriors do not agree!")
        
        # Return poisson estimates
        estimates = list(
            pi.hat = pi.hat
          , N.hat  = N.hat
          , d.hat  = d.hat
          , g.hat  = g.hat
          , X.train.parameter = s.estimate$X.train.parameter
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
        fitted.posteriors <- fit.posteriors(posteriors, levels(y))

        # Compare with Fisher's discriminant rule to assign classes
        if (type == "linear") 
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                linear.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat, pi.hat[k])))))
        if (type == "quadratic")
            discriminants <- t(apply(X, 1, function(x) unlist(lapply(1:K, function(k)
                quadratic.discriminant(as.matrix(x), mu.hat[k, ], sigma.hat[[k]], pi.hat[k])))))
        
        # Ensure we make the same predictions.
        fitted.discriminants <- fit.posteriors(discriminants, levels(y))
        
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
          , size.factor = size.factor
          , prior = prior
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
