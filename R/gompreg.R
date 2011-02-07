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
        
        R <- function(T0, T, alpha, gamma){
            ret <- exp(alpha + gamma) *
                (exp(T * exp(-gamma)) - exp(T0 * exp(-gamma)))
            ##cat("R = ", ret, "\n")
            ret
        }
        for (i in 1:ns){
            T0 <- enter[strata == i]
            T <- exit[strata == i]
            D <- event[strata == i]
            z <- X[strata == i, ,drop = FALSE]

            ezb <- exp(bz[strata == i])
            alpha <- beta[ncov + 2 * i]
            gamma <- beta[ncov + 2 * i - 1]
            driv <- D - R(T0, T, alpha, gamma) * ezb
            grad[ncov + 2 * i] <- sum(driv)
            ##cat("grad[ncov + 2 * i] = ", grad[ncov + 2 * i], "\n")
            grad[ncov + 2 * i - 1] <- -sum(D * T * exp(-gamma)) -
                sum(R(T0, T, alpha, gamma) -
                    exp(alpha) * (T * exp(T * exp(-gamma)) -
                                  T0 * exp(T0 * exp(-gamma))) * ezb)
            ##cat("grad[ncov + 2 * i - 1] = ", grad[ncov + 2 * i - 1], "\n")
            if (ncov){
                for (j in 1:ncov){
                    grad[j] <- grad[j] + sum(driv * z[, j])
                }
            }
        }
        ##cat("grad = ", grad, "\n")
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
            ##cat("shape = ", shape, ", scale = ", scale, "\n")
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
            ##cat("S1 = ", S1[1], ", S0 = ", S0[1], "\n") 
        }
        ##cat("Fmin = ", total, "\n")
        return(total)
    }

    beta <- numeric(bdim)
    if (ncov)
        beta[1:ncov] <- init  # Start values
    beta[ncov + 1] <- log(1000)
    res <- optim(beta, Fmin, gr = dGomp,
                 method = "BFGS",
                 control = list(fnscale = -1, reltol = 1e-10),
                 hessian = TRUE)
    ##cat("AFTER optim: grad = ", dGomp(res$par), "\n\n") 
    if (res$convergence != 0) stop("[gompreg]: No convergence")
    ##vari <- solve(-res$hessian)
    coefficients <- res$par
    names(coefficients) <- c(colnames(X), "log(scale)", "log(shape)")
    fit <- list(coefficients = coefficients,
                ##sd = sqrt(diag(vari)),
                loglik = res$value
                )
    fit$pfixed <- FALSE

    fit$var <- tryCatch(solve(-res$hessian), error = function(e) e)
    fit$hessian <- res$hessian
    fit$n.strata <- ns
    ##class(fit) <- "phreg"
    fit$fail <- FALSE
    fit
}
                      
