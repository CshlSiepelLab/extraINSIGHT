## Vectorized function for computing the uncertainty of a ratio
## of two estimated values X and Y given their standard errors and
## assuming cov(X,Y) = 0
ratio_estimate <- function(x, x_se, y, y_se) {
    ratio <- x / y
    ## Formula for computing SE of X/Y assuming cov(X,Y) = 0
    se <- abs(ratio) * sqrt((x_se / x)^2 + (y_se / y)^2)
    return(list(ratio, se))
}
