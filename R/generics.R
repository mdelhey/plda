plot.plda <- function(result, Xnew, y, n.seq = 100, threshold = 5, ...) {    
    N <- nrow(Xnew)
    P <- ncol(Xnew)
    is.pca <- P != 2
    stopifnot(result$parameters$P == P)

    # Use first 2 PC's to plot
    if (is.pca) {
        pca <- prcomp(Xnew, center = TRUE, scale. = TRUE)
        X.pca <- cbind(pca$x[, 1], pca$x[, 2])
        fit.pca <- plda(X.pca, y, type = result$parameters$type, prior = result$parameters$prior)
        X <- X.pca
        model <- fit.pca
    } else {
        X <- Xnew
        model <- result
    }
    
    # Generate boundries
    Xmin <- min(X[, 1]) - threshold
    Xmax <- max(X[, 1]) + threshold
    Ymin <- min(X[, 2]) - threshold
    Ymax <- max(X[, 2]) + threshold

    # Generate grid       
    X.grid <- point.grid(Xmin, Xmax, Ymin, Ymax, n.seq)
    y.grid <- predict(model, X.grid)

    # Dataframes
    df.data <- data.frame(x1 = X[, 1], x2 = X[, 2], y = y)
    df.grid <- data.frame(x1 = X.grid[, 1], x2 = X.grid[, 2], y = y.grid)

    # Plot deciscion boundry
    ggplot2::ggplot() +
      ggplot2::geom_jitter(data = df.data, ggplot2::aes(x = x1, y = x2, color = y)) +
      ggplot2::geom_point(data = df.grid, ggplot2::aes(x = x1, y = x2, color = y), alpha = 0.2) +
      ggplot2::theme_bw() +
      ggplot2::xlab("feature 1") +
      ggplot2::ylab("feature 2") + 
      ggplot2::xlim(Xmin, Xmax) +
      ggplot2::ylim(Ymin, Ymax) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
}

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

    if (result$parameters$type == "poisson") {
        densities <- estimate.poisson.densities(
            Xnew
          , lambda.hat = result$estimate$lambda.hat
        )        
    }
    
    if (result$parameters$type == "poisson.seq") {
        s.estimate <- estimate.size.factor(
            Xnew
          , X.train.parameter = result$estimates$X.train.parameter
          , type = result$parameters$size.factor
        )

        s.hat <- s.estimate$s.hat
        N.hat <- matrix(s.hat) %*% matrix(result$estimates$g.hat, nrow = 1)
        
        densities <- estimate.poisson.seq.densities(
            Xnew
          , N.hat = N.hat
          , d.hat = result$estimates$d.hat
        )
    }

    posteriors <- bayes.rule(Xnew, densities, result$estimates$pi.hat)
    fitted.posteriors <- fit.posteriors(posteriors, result$parameters$levels)

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
        printCoefmat(means, has.Pvalue = FALSE)

        cat("\nd.hat:\n")
        printCoefmat(result$estimates$d.hat, has.Pvalue = FALSE)
    }

    cat("\n")
}
