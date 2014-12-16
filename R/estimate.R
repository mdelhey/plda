#' @title Bayes rule for estimating posteriors.
#' @inheritParams plda
#' @param type Calculate posteriors using raw products of probabilities or summation of logs of
#' probabilities? "Sum" may offer superior floating point accuracy.
bayes.rule <- function(X, densities, pi, type = c("sum", "prod")) {
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
        warning("posteriors do not sum to one!")

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
    stopifnot(nrow(X) == nrow(N.hat))
    K <- nrow(d.hat)
    N <- nrow(X)
    
    densities <- do.call("rbind", lapply(1:N, function(i) unlist(lapply(1:K, function(k)
        prod(poisson.density(X[i, ], lambda = N.hat[i, ] * d.hat[k, ]))))))
    
    return(densities)
}

#' @title Estimate size factor.
#' @inheritParams plda
#' @param X.train.parameter Statistic of the training data used in calculating new s.hats. Rather
#' than calculate this from scratch using the training data, just remember it for later. This
#' parameter is different for each type.
estimate.size.factor <- function(X, X.train.parameter = NULL, type = c("mle", "quantile", "medratio")) {
    type <- match.arg(type)
    is.train <- is.null(X.train.parameter)

    x.i <- switch(
        type
      , mle = rowSums(X)
      , quantile = pmax(apply(X, 1, quantile, 0.75), 1) # Minimum quantile value is 1
      , medratio = apply(X, 1, function(x) median(x / apply(X, 2, geo.mean)))            
    )

    s.hat <- x.i / ifelse(is.train, sum(x.i), X.train.parameter)

    if (is.train && !all.equal(sum(s.hat), 1.0)) # Not true for test s.hat
        stop("Error estimating s.hat, doesn't sum to 1.")        
    
    # Need to use training data for test data.
    X.train.parameter <- sum(x.i)

    estimate <- list(
        s.hat = s.hat
      , X.train.parameter = X.train.parameter
    )

    return(estimate)
}

#' @title Fit posteriors.
fit.posteriors <- function(posteriors, levels) {
    stopifnot(is.character(levels), is.vector(levels), is.matrix(posteriors))
    factor(levels[apply(posteriors, 1, which.max)], levels = levels)
}

## #' @title Calculate discriminants
## #' @description Wrapper for discriminant functions. Not currently in use.
## calculate.discriminants <- function(X, parameters, type = c("linear", "quadratic", "poisson")) {
##     type <- match.arg(type)
##     K <- ncol(X)
##     N <- nrow(X)   
##     discriminants <- do.call("rbind", lapply(1:N, function(i) unlist(lapply(1:K, function(k)
##         poisson.discriminant(X[i, ], s.hat[i], g.hat, d.hat[k, ], pi.hat[k])))))    
##     return(discriminants)
## }
