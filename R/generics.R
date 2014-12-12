predict.plda <- function(result, Xnew, ...) {
    stopifnot(is.matrix(Xnew), ncol(Xnew) == result$parameters$P)
    stopifnot(class(result) == "plda")

    if (result$parameters$type %in% c("linear", "quadratic")) {
        
        densities <- estimate.normal.densities(
            Xnew
          , mu.hat = result$estimates$mu.hat
          , sigma.hat = result$estimates$sigma.hat
          , type = result$parameters$type
        )
    }

    posteriors <- bayes.rule(Xnew, densities, result$estimates$pi.hat)

    fitted.posteriors <- factor(result$parameters$levels[apply(posteriors, 1, which.max)],
                                levels = result$parameters$levels)

    return(fitted.posteriors)
}

print.plda <- function(result, ...) {
    dist <- ifelse(result$parameters$type == "poisson", "poisson", "normal")
    cat(sprintf("X(%i x %i) with K = %i\n", result$parameters$N, result$parameters$P, result$parameters$K))
    cat(sprintf("Type: **%s** discriminant analysis", result$parameters$type), "\n")
    cat(sprintf("Dist: predictors have a **%s** distribution", dist), "\n")

    cat("\npi.hat:\n")
    printCoefmat(t(as.matrix(result$estimates$pi.hat)))

    if (dist == "normal") {
        cat("\nmu.hat:\n")
        printCoefmat(result$estimates$mu.hat)

        cat("\nsigma.hat:\n")
        printCoefmat(result$estimates$sigma.hat)
    }

    if (dist == "poisson") {
        cat("\nN.hat (N by p, summarized by column averages):\n")
        means <- t(as.matrix(colMeans(as.matrix(result$estimates$N.hat))))
        printCoefmat(means)

        cat("\nd.hat:\n")
        printCoefmat(result$estimates$d.hat)
    }

    cat("\n")
}

