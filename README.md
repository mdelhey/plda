plda
====

poisson discriminant analysis 

includes routines for gaussian LDA and QDA 

## install

```{r}
install.packages("devtools")
devtools::install_github("mattdelhey/plda")
```

## useage

```{r}
library(plda)
data(iris)
X <- as.matrix(iris[, 1:4])
y <- as.factor(iris[, 5])
fit <- plda(X, y, method = "linear", type = "normal")
fit
fitted(fit)
```

