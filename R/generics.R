predict.plda <- function(result, Xnew, ...) {
    stopifnot(is.matrix(Xnew), ncol(Xnew) == result$parameters$P)   
    #estimate.densities(result)
}

print.plda <- function(result, ...) {
    cat(sprintf("X(%i x %i) with K = %i\n", result$parameters$N, result$parameters$P, result$parameters$K))
    cat("Method: %s discriminant analysis\n", result$parameters$method)
    cat("Dist: predictors have %s distribution\n", result$parameters$dist)

    cat("\npi.hat:\n")
    printCoefmat(t(as.matrix(result$pi.hat)))

    cat("\nmu.hat:\n")
    printCoefmat(result$mu.hat)    
}

