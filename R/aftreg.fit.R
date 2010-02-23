aftreg.fit <- function(X, Y, dist,
                       strata, offset,
                       init, shape, id,
                       control, center = FALSE){

    ## New in Version 1.2-9; wrong before!
    ## Note that we MUST keep individuals together here;
    ## stratum comes second, beacuse we can then let individuals
    ## change stratum over time!
    ord <- order(id, Y[, 1])
    X <- X[ord, , drop = FALSE]
    Y <- Y[ord, , drop = FALSE]
    id <- id[ord]
    strata <- strata[ord]
    offset <- offset[ord]
    #####################################
    if (dist == "weibull"){
        dis <- 0
    }else if(dist == "loglogistic"){
        dis <- 1
    }else if (dist == "lognormal"){
        dis <- 2
    }else if (dist == "ev"){
        dis <- 3
    }else if (dist == "gompertz"){ # An EV with shape == 1:
        ## dis <- 4
        dis <- 3
        shape <- 1
    }else{
        stop(paste(dist, "is not an implemented distribution"))
    }

    nn <- NROW(X)
    ncov <- NCOL(X)

    ## No intercepts in 'aftreg!!! (1.2-17) intercept <- (dis == 4) # gompertz
    if (ncov){
        wts <- Y[, 2] - Y[, 1]
        means <- apply(X, 2, weighted.mean, w = wts)
    ##    if (intercept) means[1] <- 0
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


    ## Not needed?? nstra <- c(0, cumsum(table(strata)))
    if (all(shape <= 0)){ ## Then shape is estimated in all strata

        Fmin <- function(beta){

            fit <- .C("aftsup",
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
        ## Start values (no covariates):
        bdim <- 2 * ns
        ncov.save <- ncov
        ncov <- 0
        beta <- numeric(bdim)
        for (i in 1:ns){
            beta[2 * i - 1] <- log(sum(Y[, 2] - Y[, 1]) / sum(Y[, 3]))
            beta[2*i] <- 0
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
                    pi.hat <- exp(fit$beta[col])
                    scale.corr <- sum(means * fit$beta[1:ncov]) /
                        pi.hat
                    fit$beta[row] <- fit$beta[row] + scale.corr
                    dxy[row, 1:ncov] <- means / pi.hat
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


        fit$shape.fixed <- FALSE

    }else{  ## Then shape is fixed in all strata:
        ## Note: We must allow stratification even here (091006)!!


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
        if (ncov & !center){
            dxy <- diag(bdim)
            dxy[bdim, 1:ncov] <- means / shape
            scale.corr <- sum(means * fit$beta[1:ncov]) / shape
            fit$beta[bdim] <- fit$beta[bdim] + scale.corr
            dxy[bdim, bdim] <- -1
        }
        coef.names <- c(colnames(X), "log(scale)")

        ## Note; this is really a "hack"!!!!!!!!!!!!!!!
    }

    ##cat("done!\n")

    if (!fit$fail){
        if (ncov && is.numeric(fit$variance)){
            var <- dxy %*% matrix(fit$variance, bdim, bdim) %*% t(dxy)
            colnames(var) <- rownames(var) <- coef.names
        }else{
            var <- fit$variance
        }
    }
    else
        var <- NULL

    coefficients <- fit$beta
    names(coefficients) <- coef.names

    list(coefficients = coefficients,
         df = ncov,
         var = var,
         loglik = fit$loglik,
         score = fit$sctest,
         conver = fit$conver,
         fail = fit$fail,
         iter = fit$iter,
         n.strata = ns,
         shape = fit$shape
         )

}
