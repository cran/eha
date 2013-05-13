gompstartAft <- function(enter, exit, event){
    ## Gives start values for Gompertz parameters
    ## To be used in aftreg.
    ##
    ## This is for ONE stratum only! So input is from only one!

    ## Profiling; gamma is profiled out:
    D <- sum(event)
    logD <- log(D)
    DT <- sum(exit * event)
    funk <- function(alpha){
        n <- length(alpha) # Vectorizing
        loglik <- numeric(n)
        for (i in 1:n){
            ealpha <- exp(-alpha[i])
            S <- sum(exp(exit * ealpha) - exp(enter * ealpha))
            gamma <- logD - alpha[i] - log(S)
            loglik[i] <- D * gamma + DT * ealpha - D
        }
        loglik
    }

    alpha <- log(max(exit)) # start value
    ##from <- alpha - width / 2
    ##to <- alpha + width / 2
    fit <- optim(alpha, funk, control = list(fnscale = -1), method = "BFGS")
    alpha <- fit$par
    S <- sum(exp(exit * exp(-alpha)) - exp(enter * exp(-alpha)))
    gamma <- log(D) - alpha - log(S)
    ret <- c(alpha, gamma)
    ##names(ret) <- c("alpha", "gamma")
    ret
}
    
