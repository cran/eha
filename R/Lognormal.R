hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), log = FALSE){
    ## shape = 1 / sdlog, scale = exp(meanlog)
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- dlnorm(x, meanlog, sdlog) /
        plnorm(x, meanlog, sdlog, lower.tail = FALSE)
    if (log) ret <- ifelse(ret <= 0, -Inf, log(ret))
    return (ret)

}

Hlnorm <- function(x, meanlog = 0, sdlog = 1,
                   shape = 1 / sdlog, scale = exp(meanlog), log.p = FALSE){
    meanlog <- log(scale)
    sdlog <- 1 / shape
    ret <- -plnorm(x, meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
    if (log.p) ret <- log(ret)
    return (ret)
}
