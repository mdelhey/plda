predict.plda <- function(result, Xnew, ...) {
    stopifnot(is.matrix(Xnew), ncol(Xnew) == result$parameters$P)
    stopifnot(class(result) == "plda")
    
    densities <- estimate.densities(
        Xnew
      , mu.hat = result$estimates$mu.hat
      , sigma.hat = result$estimates$sigma.hat
      , method = result$parameters$method
    )

    posteriors <- bayes.rule(Xnew, densities, result$estimates$pi.hat)

    fitted.posteriors <- factor(result$parameters$levels[apply(posteriors, 1, which.max)],
                                levels = result$parameters$levels)

    return(fitted.posteriors)
}

print.plda <- function(result, ...) {
    cat(sprintf("X(%i x %i) with K = %i\n", result$parameters$N, result$parameters$P, result$parameters$K))
    cat(sprintf("Method: **%s** discriminant analysis", result$parameters$method), "\n")
    cat(sprintf("Dist: predictors have **%s** distribution", result$parameters$dist), "\n")

    cat("\npi.hat:\n")
    printCoefmat(t(as.matrix(result$pi.hat)))

    cat("\nmu.hat:\n")
    printCoefmat(result$mu.hat)

    cat("\n")
}

