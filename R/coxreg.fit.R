coxreg.fit <- function(X, Y, rs, strats, offset, init, max.survs,
                       method = "breslow", boot = FALSE, efrac = 0,
                       calc.hazards = TRUE, calc.martres = TRUE,
                       control, verbose = TRUE){

    nn <- NROW(X)
    ncov <- NCOL(X)

    if (missing(strats) || is.null(strats)) 
      strats <- rep(1, nn)

    if (missing(rs) || is.null(rs)){
        rs <- risksets(Y, strata = strats, max.survs)
    }

    
    if (max(rs$riskset) > nn) stop("Riskset does not match data")
    
    if (missing(offset) || is.null(offset)) 
      offset <- rep(0, nn)
    
    if (missing(init) || is.null(init)) 
      init <- rep(0, ncov)
    
    if (missing(control)){
        control <- list(eps=1.e-8, maxiter = 10, trace = FALSE)
    }else{
        if (!is.numeric(control$eps)){
            stop("Error in control = list(eps = ...) ")
        }else{
            if (control$eps <= 0) stop("control$eps must be strictly positive")
        }
        if (!is.numeric(control$maxiter)){
            stop("Error in control = list(maxiter = ...) ")
        }else{
            if (control$maxiter < 0) stop("control$maxiter must be positive")
        }
        if (!is.logical(control$trace)) stop("control$trace must be logical")
    }
    
    nullModel <- NCOL(X) == 0
    
    if (nullModel){
        cat("nullModel!!!!!!!!!!!!!!!!!!!!!\n")
        ## faking a simple model with no iterations
        ncov <- 0
        X <- matrix(0, ncol = 1, nrow = nn)
        control$maxiter <- 0
        init <- 0
    }
    
    printlevel <- control$trace
    ## NOTE: silent == TRUE ===> printlevel = 0
    iter <- control$maxiter
    if (method[1] == "efron")
      meth <- 0
    else if (method[1] == "breslow")
      meth <- 1
    else if (method[1] == "mppl")
      meth <- 2
    else if (method[1] == "ml")
      meth <- 3
    else
      stop(paste("Unknown method", as.character(method[1]))) 

    boot <- abs(as.integer(boot))

    if (!nullModel){
        fit <- .C("sup",
                  as.integer(meth),
                  iter = as.integer(iter), #maxit on input, actual on output
                  as.double(control$eps),
                  as.integer(printlevel),
                                        #
                  as.integer(sum(rs$n.events)), ## total No. of events
                  as.integer(sum(rs$antrs)),  ## total No. of risksets
                  as.integer(length(rs$antrs)), # No. of strata
                                        #
                  as.integer(rs$antrs),
                  as.integer(rs$n.events),
                  as.integer(rs$size),
                                        #
                  as.integer(length(rs$riskset)), # Sum of risk set sizes.
                  as.integer(rs$eventset),
                  as.integer(rs$riskset),
                                        #
                  as.integer(nn),
                  as.integer(ncov),
                  ## Note here; X is transposed! From 2007-04-16 (0.99)
                  as.double(t(scale(X, center = TRUE, scale = FALSE))),
                  as.double(offset),
                                        #
                  as.double(init),     # 'start.beta'
                  boot = as.integer(boot),
                  as.double(efrac),
                  beta = double(ncov * (1 + boot)),
                  sd.beta = double(ncov * (1 + boot)),
                                        #
                  loglik = double(2), # [1] == start, [2] == maximized
                  variance = double(ncov * ncov),
                  sctest = double(1),
                  ##
                  hazard = double(sum(rs$antrs)),
                  conver = integer(1),
                  f.conver = integer(1),
                  fail = integer(1),
                  DUP = FALSE,
                  PACKAGE = "eha")
    
        bootstrap <- NULL
        boot.sd <- NULL
        if (boot & (fit$fail == 0)){
            
            bootstrap <- matrix(fit$beta[(ncov + 1):((boot + 1) * ncov)],
                                nrow = ncov, ncol = boot)
            boot.sd <- matrix(fit$sd.beta[(ncov + 1):((boot + 1) * ncov)],
                              nrow = ncov, ncol = boot)
            fit$beta <- fit$beta[1:ncov]
            fit$sd.beta <- fit$sd.beta[1:ncov]
        }   
    }else{ # if nullModel
        X <- matrix(0, ncol = 0, nrow = nn)
        fit <- list(beta = numeric(0),
                    conver = TRUE,
                    fail = FALSE)
        ncov <- 0
        bootstrap <- NULL
        boot.sd <- NULL
        
    }
    if (fit$fail){
        out <- paste("Singular hessian; suspicious variable No. ",
                     as.character(fit$fail), ":\n",
                     colnames(X)[fit$fail], sep = "")
        if (verbose) warning(out)## New
        return(fit)## 19 February 2007.
    }else if (!fit$conver){
        ##fit$conver <- 1 Removed 19 February 2007
        if (!fit$f.conver){
            if (verbose)
              warning("Did not converge")
        }else{
            if (verbose)
              warning("log likelihood converged, but not variables")
        }
    }
    
    if (calc.hazards && (!fit$fail)){
        score <- exp(X %*% fit$beta)
        hazard <- .Fortran("hazards",
                           as.integer(sum(rs$antrs)),  ## total No. of risksets
                           as.integer(length(rs$antrs)), # No. of strata
                           ##
                           as.integer(rs$antrs),
                           as.integer(rs$n.events),
                           as.integer(rs$size),
                           ##
                           as.integer(length(rs$riskset)),
                           ## is Sum of risk set sizes.
                           as.integer(rs$riskset),
                           ##
                           as.integer(nn),
                           ##
                           as.double(score),
                           hazard = double(sum(rs$antrs)),
                           ##
                           DUP = FALSE,
                           PACKAGE = "eha"
                           )$hazard
    }else{ #if not calc.hazards or fail
        hazard <- NULL
    }
    
    if (calc.martres && (!nullModel)){
        resid <- .Fortran("martres",
                          as.integer(sum(rs$antrs)),
                          as.integer(length(rs$antrs)),
                          ##
                          as.integer(rs$antrs),
                          as.integer(rs$n.events),
                          as.integer(rs$size),
                          ##
                          as.integer(length(rs$riskset)), # Sum of risk set sizes.
                          as.integer(rs$riskset),
                          ##
                          as.integer(nn),
                          ##
                          as.double(score),  # 'score'
                          as.double(hazard),
                          resid = double(nn),
                          DUP = FALSE,
                          PACKAGE = "eha"
                          )$resid
    }else{
        resid <- NULL
    }
    
    if ((!fit$fail) && (!nullModel))
      var <- matrix(fit$variance, ncov, ncov)
    else
      var <- NULL


    list(coefficients = fit$beta,
         sd = fit$sd.beta,
         var = var,
         loglik = fit$loglik,
         score = fit$sctest,
         linear.predictors = X %*% fit$beta,
         residuals = resid,
         noOfRisksets = rs$antrs,
         risktimes = rs$risktimes,
         hazard = hazard,
         means = apply(X, 2, mean),
         bootstrap = bootstrap,
         boot.sd = boot.sd,
         conver = fit$conver,
         f.conver = fit$f.conver,
         fail = fit$fail,
         iter = fit$iter
         )    
}
