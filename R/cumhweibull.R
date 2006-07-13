Hweibull <- function(x, shape, scale = 1, log = FALSE){
    res <- (x / scale)^shape
    if (log) res <- logb(res, base = exp(1))
    res
}
