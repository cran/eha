hweibull <- function(x, shape, scale = 1, log = FALSE){
    if (shape <= 0 || scale <= 0)
      stop("scale and shape must be positive")  
    res <- ifelse(x < 0, 0, shape * (x / scale)^(shape - 1) / scale) 
    if (log) res <- log(x)
    res
}
