aftp0g <- function(printlevel, ns, nn, id,
                   strata, Y, X, offset, means){
    ## Gompertz aft-reg!
    Fmin <- function(beta){

        fit <- .C("aftregGomp",
                  as.integer(printlevel),
                  ##
                  as.integer(ns), # No. of strata
                  as.integer(nn),
                  as.integer(ncov),
                  as.integer(bdim),
                  ##
                  as.integer(id),
                  as.integer(strata - 1), # 0, ..., (ns-1); C-style!
                  as.double(Y[, 1]),  ## 'enter'
                  as.double(Y[, 2]),  ## 'exit'
                  as.integer(Y[, 3]), ## 'event'
                  ##
                  ##as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(t(X)), # NOTE; scaling already done!
                  ## NOTE transpose!
                  as.double(offset),
                  as.integer(dis),     # baseline distribution
                  as.double(beta),
                                        # results -->
                  loglik = double(1),  # function value at beta
                  fail = integer(1), # = 0: No failure
                  DUP = FALSE,
                  PACKAGE = "eha"
                  )
        if (fit$fail) stop("Error in likelihood calculation")
        return(fit$loglik)
    }
    ncov <- NCOL(X)
    dis <- 4 ## Just a joke to make Gomp happy
    ## Start values (no covariates):
    bdim <- 2 * ns
    ncov.save <- ncov
    ncov <- 0
    beta <- numeric(bdim)
    for (i in 1:ns){
        enter <- Y[strata == 1, 1]
        exit <-  Y[strata == 1, 2]
        event <-  Y[strata == 1, 3]
        
        alpha <- log(max(exit))
        beta[(2 * i - 1):(2 * i)] <- gompstart(enter, exit, event, width = 20)
        ##beta[2 * i - 1] <- log(sum(Y[, 2] - Y[, 1]) / sum(Y[, 3]))
        ##beta[2*i] <- 0
    }
    
    loglik.start <- -optim(beta, Fmin, method = "BFGS", hessian = FALSE)$value
    ## And now the real thing:
    ncov <- ncov.save
    bdim <- ncov + 2 * ns
    beta <- c(rep(0, ncov), beta)
    fit <- optim(beta, Fmin, method = "BFGS",
                 control = list(trace = as.integer(printlevel)),
                 hessian = TRUE)
    fit$beta <- fit$par
    fit$loglik <- c(loglik.start, -fit$value)
    fit$variance <- try(solve(fit$hessian))
    fit$fail <- FALSE
    if (ncov){
        dxy <- diag(2 * ns + ncov)
        for (i in 1:ns){ ## Really a HACK ??!!!!!!!!!!!!!!!
            row <- ncov + 2 * i - 1
            col <- row + 1
            ## fit$beta[row] <- -fit$beta[row] NOT ANY MORE!
            if (ncov){
                ##pi.hat <- exp(fit$beta[col])
                scale.corr <- sum(means * fit$beta[1:ncov]) #/
                ##pi.hat
                fit$beta[row] <- fit$beta[row] + scale.corr
                dxy[row, 1:ncov] <- means # / pi.hat
                dxy[row, col] <- -scale.corr
            }
            ##dxy[row, row] <- -1
        }
    }
    
    coef.names <- colnames(X)
    if (ns > 1){
        for (i in 1:ns){
            coef.names <- c(coef.names,
                            paste("log(scale)", as.character(i), sep =":"),
                            paste("log(shape)", as.character(i), sep =":"))
        }
        
    }else{
        coef.names <- c(coef.names,
                        "log(scale)", "log(shape)")
    }
    
    fit$coef.names <- coef.names
    fit$shape.fixed <- FALSE

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
