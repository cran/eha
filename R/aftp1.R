aftp1 <- function(printlevel, ns, nn, id,
                  strata, Y, X, offset, shape, dis, means){
    Fexpmin <- function(beta){

        fit <- .C("aftexpsup",
                  as.integer(printlevel),
                                        #
                  as.integer(ns), # No. of strata
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                                        #
                  as.integer(id),
                  as.integer(strata - 1), # 0, ..., (ns-1); C-style!
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                                        #
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t((X))), #NOTE: scaling already done!
                  as.double(offset),
                  as.double(shape), ## "p" (fixed)
                  as.integer(dis),
                  beta = as.double(beta),
                                        # results -->
                  loglik = double(1), # Return value at beta
                  fail = integer(1),
                  DUP = TRUE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) stop("Error in exp likelihood calculation")
        
        return(fit$loglik)
    }
    ncov <- NCOL(X)
    if (length(shape) == 1){
        shape <- rep(shape, ns)
    }else if (length(shape) != ns){
        stop("length(shape) must be equal to 1 or No. of strata")
    }
    
    bdim <- ns
    ncov.save <- ncov
    ncov <- 0
    beta <- numeric(bdim)
    for (i in 1:ns){
        beta[i] <- log(sum(Y[, 2] - Y[, 1]) / sum(Y[, 3]))
    }
    
    res <- optim(beta, Fexpmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = FALSE)
    ncov <- ncov.save
    bdim <- ncov + ns
    beta <- c(rep(0, ncov), res$par)
    loglik.start <- -res$value
    
    fit <- optim(beta, Fexpmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = TRUE)
    fit$fail <- (fit$convergence > 0.5)
    fit$beta <- fit$par
    fit$loglik <- c(loglik.start, -fit$value)
    fit$variance <- try(solve(fit$hessian))
    fit$shape.fixed <- TRUE
    fit$shape <- shape
    fit$shape.sd <- NULL  ## Not necessary!?!?
    ##fit$beta[bdim] <- -fit$beta[bdim] # To get "1 / lambda"! NO!!
    ##if (ncov & !center){
    if (ncov){
        dxy <- diag(bdim)
        dxy[bdim, 1:ncov] <- means# / shape
        scale.corr <- sum(means * fit$beta[1:ncov])# / shape
        fit$beta[bdim] <- fit$beta[bdim] + scale.corr
        dxy[bdim, bdim] <- -1
    }
    fit$coef.names <- c(colnames(X), "log(scale)")

    if (!fit$fail){
        if (ncov && is.numeric(fit$variance)){
            fit$var <- dxy %*% matrix(fit$variance, bdim, bdim) %*% t(dxy)
            colnames(fit$var) <- rownames(fit$var) <- fit$coef.names
        }else{
            fit$var <- fit$variance
        }
    }
    else
        fit$var <- NULL

    fit$ncov <- ncov
    fit
}
