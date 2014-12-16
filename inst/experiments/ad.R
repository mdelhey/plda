library(plda)
library(ggplot2)
set.seed(1)

# For development
devtools::install("~/plda")
devtools::load_all("~/plda")

################################################################################
### Ad data
################################################################################

fname <- "~/plda/inst/extdata/ad/ad.data"
ad <- read.csv(fname, header = FALSE, stringsAsFactors = FALSE)

ad[which(stringr::str_detect(ad[, 1], "?")), 1]
ad[which(stringr::str_detect(ad[, 1], "?")), 1]
x <- stringr::str_detect(ad[, 1], "?")

ad[]
y <- as.factor(data[, ncol(data)])
X <- as.data.frame(data[, -ncol(data)])
X[, 1:3] == ?
X[1:100, 1:3]
X[73,1]

