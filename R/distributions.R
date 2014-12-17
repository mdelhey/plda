#' @title Multivariate normal density.
normal.density <- function(x, mu, sigma) {
    det(2*pi*sigma)^(-1/2) * exp(-1/2 * t(x - mu) %*% solve(sigma) %*% (x - mu))
}

#' @title Poisson univariate density.
poisson.density <- function(x, lambda) {
    (lambda^x / factorial(x)) * exp(1)^(-lambda)
}

#' @title Linear discriminant.
linear.discriminant <- function(x, mu, sigma, pi) {    
    t(x) %*% solve(sigma) %*% mu - 0.5*t(mu) %*% solve(sigma) %*% mu + log(pi)
}

#' @title Quadratic discriminant.
quadratic.discriminant <- function(x, mu, sigma, pi) {
    (-1/2)*log(det(sigma)) - (1/2)*t(x - mu) %*% solve(sigma) %*% (x - mu) + log(pi)
}

#' @title Poisson sequence discriminant.
poisson.seq.discriminant <- function(x, s, g, d, pi) {
    sum(x*log(d)) - s*sum(g*d) + log(pi)
}

#' @title Poission discriminant.
poisson.discriminant <- function(x, lambda, pi) {
    sum(x*log(lambda) - log(factorial(x)) - lambda) + log(pi)
}
