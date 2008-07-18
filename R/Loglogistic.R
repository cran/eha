## Contains functions for the log-logistic distribution

pllogis <- function(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    location <- log(scale)
    scale <- 1 / shape
    if (log.p){
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          -Inf,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          0,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }
    }else{ #! log.p
        if (lower.tail){
            ret <- ifelse(q <= 0,
                          0,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }else{
            ret <- ifelse(q <= 0,
                          1,
                          plogis(log(q), location, scale, lower.tail, log.p)
                          )
        }
    }

    return ( ret )
}

dllogis <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    y <- x / scale


    ret <- ifelse(y < 0,
                  0,
                  shape * y^(shape - 1) / (1 + y^shape)^2
                  )
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))

    return (ret)
}

hllogis <- function(x, shape = 1, scale = 1, log = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    y <- x / scale
    ret <- ifelse(x < 0,
                    0,
                    shape * y^(shape - 1) / (1 + y^shape)
                    )
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))

    return(ret)
}

qllogis <- function(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    if (log.p) p <- exp(p)

    location <- log(scale)
    scale <- 1 / shape

    ## Let 'qlogis do the checking....
    ret <- qlogis(p, location, scale, lower.tail, log.p)

    return(exp(ret))
}



Hllogis <- function(x, shape = 1, scale = 1, log.p = FALSE){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }

    ret = -pllogis(x, shape, scale, lower.tail = FALSE, log.p = TRUE)
    if (log.p) ret <- log(ret)

    return (ret)
}

rllogis <- function(n, shape = 1, scale = 1){
    if ( any(c(shape, scale) <= 0) ){
        warning("Non-positive shape or scale")
        return(NaN)
    }
    location <- log(scale)
    scale <- 1 / shape

    exp(rlogis(n, location, scale))
}
