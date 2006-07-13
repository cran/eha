hweibull <- function(x, shape, scale = 1, log = FALSE){
    if (shape <= 0 || scale <= 0)
      error("scale and shape must be positive")  
    (x / scale)^(shape - 1) / scale
}
