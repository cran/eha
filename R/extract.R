extractAIC.coxreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

extractAIC.phreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

extractAIC.aftreg <- function(fit, scale, k = 2, ...) 
{
    edf <- sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}

nobs.coxreg <- function(object, ...){
    object$n
}

nobs.phreg <- function(object, ...){
    object$n
}
nobs.aftreg <- function(object, ...){
    object$n
}
