plda
====

Poisson discriminant analysis with additional routines for gaussian LDA and QDA.

See [Classification and lustering of sequencing data using a poisson model](http://arxiv.org/pdf/1202.6201.pdf).

## Installation

```{r}
install.packages("devtools")
devtools::install_github("mattdelhey/plda")
```

## Useage

```{r}
library(plda)
data(iris)
X <- as.matrix(iris[, 1:4])
y <- as.factor(iris[, 5])
fit <- plda(X, y, type = "quadratic")
fit
fitted(fit)
```
