gompreg <- function(X, Y, strata, offset, init, control){
    ## Gompertz proportional hazards:
    ## Stratum i: h_i(t; beta) = a_i * exp(t / b_i) * exp(x*beta),
    ## i = 1, ..., ns
    ## NOTE: Y is nn x 3!
    
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    nn <- NROW(X)
    ncov <- NCOL(X)
    if (NROW(Y) != nn) stop("Y NROW error")
    if (NCOL(Y) != 3) stop("Y NCOL error", )
    if (ncov){
        wts <- Y[, 2] - Y[, 1]
        means <- apply(X, 2, weighted.mean, w = wts)
        means <- apply(X, 2, mean)

        for (i in 1:ncov){
            X[, i] <- X[, i] - means[i]
        }
    }
    if (missing(strata) || is.null(strata)){
        strata <- rep(1, nn)
        ns <- 1
    }else{
        strata <- as.integer(factor(strata))
        ns <- max(strata)
    }

    if (length(strata) != nn) stop("Error in stratum variable")
    if (missing(offset) || is.null(offset))
        offset <- rep(0, nn)

    if (missing(init) || is.null(init))
        init <- rep(0, ncov)
    if (length(init) != ncov) stop("Error in init")

    printlevel <- control$trace
    iter <- control$maxiter


    nstra <- c(0, cumsum(table(strata)))

    bdim <- ncov + 2 * ns

    dGomp <- function(beta){
        ## Calculates the first derivatives of a Gompertz regression (stratified)
        ## beta[1:ncov] = the beta coefficients.
        ## beta[ncov + 1], beta[ncov + 3], ... = scale[1, 2, ...] ("gamma")
        ## beta[ncov + 2], beta[ncov + 4], ... = shape[1, 2, ...] ("alpha")


        enter <- Y[, 1]
        exit <- Y[, 2]
        event <- Y[, 3]
        if (ncov){
            b <- beta[1:ncov]
            bz <- offset + X %*% b
        }else{
            bz <- offset
        }
        
        grad <- numeric(ncov + 2 * ns)
        
        ##RR <- function(T0, T, alpha, gamma){ # Changed 1.3-1; moved gamma;
          ##  ret <- exp(alpha) *
            ##    
        ##    ret
        ##}
        for (i in 1:ns){
            T0 <- enter[strata == i]
            T <- exit[strata == i]
            D <- event[strata == i]
            z <- X[strata == i, ,drop = FALSE]

            ezb <- exp(bz[strata == i])
            alpha <- beta[ncov + 2 * i]
            eapzb <- exp(alpha + bz[strata == i])
            gamma <- beta[ncov + 2 * i - 1]
            emgamma <- exp(-gamma)
            R <- exp(gamma + T * emgamma) - exp(gamma + T0 * emgamma) 
            driv <- D - R * eapzb
            grad[ncov + 2 * i] <- sum(driv)
            grad[ncov + 2 * i - 1] <- -sum(D * T) * emgamma -
                sum(eapzb * (R -
                    (T * exp(T * emgamma) - T0 * exp(T0 * emgamma))))
            if (ncov){
                for (j in 1:ncov){
                    grad[j] <- grad[j] + sum(driv * z[, j])
                }
            }
        }
        grad
    }    
    

    Fmin <- function(beta){
        total <- 0
        if (ncov){
            b <- beta[1:ncov]
        }else{
            b <- 0
        }
        for (i in 1:ns){
            scale <- exp(beta[ncov + 2 * i - 1])
            shape <- exp(beta[ncov + 2 * i])
            if (ncov){
                bz <- offset[strata == i] + X[strata == i, , drop = FALSE] %*% b
            }else{
                bz <- offset[strata == i]
            }
            ebz <- exp(bz)
            S1 <- pgompertz(Y[strata == i, 2], shape = shape, scale = scale,
                        log.p = TRUE, lower.tail = FALSE)
            S0 <- pgompertz(Y[strata == i, 1], shape = shape, scale = scale,
                        log.p = TRUE, lower.tail = FALSE)
            h <- hgompertz(Y[strata == i, 2], shape = shape, scale = scale,
                           log = TRUE)
            ret1 <- sum(Y[strata == i, 3] * (h + bz))
            ret2 <- sum(ebz * (S1 - S0))
            total <- total + ret1 + ret2
        }

        return(total)
    }

    ## First, fit the 'null' model:
    beta0 <- numeric(2 * ns)
    ncov.save <- ncov
    ncov <- 0
    ## Start values (primitive!!): # Better now (1.4)
    for (i in 1:ns){
        enter <- Y[strata == i, 1]
        exit <- Y[strata == i, 2]
        event <- Y[strata == i, 3]
        beta0[(ncov + 2 * i - 1):(ncov + 2 * i)] <-
            gompstart(enter, exit, event, width = 20) 
        ##beta0[ncov + 2 * i - 1] <- log(max(Y[strata == i, 2]))
        ##beta0[ncov + 2 * i] <- log(sum(Y[strata == i, 3]) /
          ##                          sum(Y[strata == i, 2])) - 1
    }

    res0 <- optim(beta0, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10),
                 hessian = FALSE)
    ## Done; now the real thing:
    ncov <- ncov.save
    beta <- numeric(bdim)
    
    if (ncov)
        beta[1:ncov] <- init  # Start values
    beta[(ncov + 1):length(beta)] <- res0$par # Ditto
    res <- optim(beta, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10),
                 hessian = TRUE)
    if (res$convergence != 0) stop("[gompreg]: No convergence")
    coefficients <- res$par
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

    names(coefficients) <- coef.names
    fit <- list(coefficients = coefficients,
                loglik = c(res0$value, res$value)
                )
    fit$pfixed <- FALSE
    fit$var <- tryCatch(solve(-res$hessian), error = function(e) e)
    if (is.matrix(fit$var)){
        colnames(fit$var) <- coef.names
        rownames(fit$var) <- coef.names
    }
    fit$hessian <- res$hessian
    if (is.matrix(fit$hessian)){
        colnames(fit$hessian) <- coef.names
        rownames(fit$hessian) <- coef.names
    }
    
    fit$n.strata <- ns
    fit$df <- ncov
    fit$fail <- FALSE # Optimist!
    fit
}
                      
