### Contains functions for the Gompertz distribution

pgompertz <- function(q, shape = 1, scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ret <- pEV(q, 1, scale, lower.tail, log.p)^shape

    return ( ret )
}

dgompertz <- function(x, shape = 1, scale = 1, log = FALSE){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ret <- dEV(x, 1, scale, log = FALSE) * shape *
        pEV(x, 1, scale,
                  lower.tail = TRUE, log.p = FALSE)^(shape - 1)
    if (log) ret <- log(ret)

    return ( ret )
}

hgompertz <- function(x, shape = 1, scale = 1, log = FALSE){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    ret <- shape * hEV(x, 1, scale, log = FALSE)
    if (log) ret <- log(ret)

    return ( ret )
}

qgompertz <- function(p, shape = 1, scale = 1,
                     lower.tail = TRUE, log.p = FALSE){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    if (log.p) p <- exp(p)

    ok <- (p >= 0) & (p <= 1)


    ret <- ifelse(ok,
                  scale * ( log1p(-log(p) / shape) ),
                  NaN)
    if (!all(ok)) warning("qgompertz produced NaN's")

    return ( ret )
}

Hgompertz <- function(x, shape = 1, scale = 1, log.p = FALSE){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    ret <- HEV(x, 1, scale, log.p = FALSE) * shape
    if (log.p) ret <- log(ret)

    return ( ret )
}

rgompertz <- function(n, shape = 1, scale = 1){

    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    y <- runif(n)

    return ( qgompertz(y, shape, scale) )
}
