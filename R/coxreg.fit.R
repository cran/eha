# Cox regression. (C) Göran Broström (2003). Initial code borrowed from
# 'coxph' in 'survival. Thanks to Terry Therneau and Thomas Lumley. 

coxreg.fit <- function(X, Y, rs, strats, offset, init, max.survs,
                       method = "breslow", boot = FALSE, control){

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
      

  printlevel <- control$trace
      ## NOTE: silent == TRUE ===> printlevel = 0
  iter <- control$maxiter
  fit <- .Fortran("sup",
                  as.integer(method == "efron"),
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
                  as.double(scale(X, center = TRUE, scale = FALSE)),
                  as.double(offset),
                  #
                  as.double(init),     # 'start.beta'
                  beta = double(ncov),
                  #
                  loglik = double(2), # [1] == start, [2] == maximized
                  dloglik = double(ncov),
                  variance = double(ncov * ncov),
                  sctest = double(1),
                  #
                  double(nn),     ## 'score', work area
                  double(ncov),          ## 'sumdscore', work area.
                  double(ncov * ncov), ## 'sumd2score', work area.
                  #
                  conver = integer(1),
                  f.conver = integer(1),
                  fail = integer(1),
                  DUP = FALSE)

  if (fit$fail){
      out <- paste("Singular hessian; suspicious variable No. ",
                   as.character(fit$fail), ":\n",
                   colnames(X)[fit$fail], sep = "")
      warning(out)## New
      return(NULL)## 5 May 2004.
  }else if (!fit$conver){
      fit$conver <- 1
      if (!fit$f.conver){
          warning("Did not converge")
      }else{
          warning("log liklihood converged, but not variables")
      }
  }
      
  score <- exp(X %*% fit$beta)
  hazard <- .Fortran("hazards",
                     as.integer(sum(rs$n.events)), ## total No. of events
                     as.integer(sum(rs$antrs)),  ## total No. of risksets
                     as.integer(length(rs$antrs)), # No. of strata
                                        #
                     as.integer(rs$antrs),
                     as.integer(rs$n.events),
                     as.integer(rs$size),
                                        #
                     as.integer(length(rs$riskset)), # Sum of risk set sizes.
                     as.integer(rs$riskset),
                                        #
                     as.integer(nn),
                     as.integer(ncov),
                                        #
                     as.double(fit$beta),
                     as.double(score),
                     hazard = double(sum(rs$antrs))
                                        #
                     #DUP = FALSE
                     )$hazard

  resid <- .Fortran("martres",
                    as.integer(sum(rs$n.events)),
                    as.integer(sum(rs$antrs)),
                    as.integer(length(rs$antrs)),
                    #
                    as.integer(rs$antrs),
                    as.integer(rs$n.events),
                    as.integer(rs$size),
                    #
                    as.integer(length(rs$riskset)), # Sum of risk set sizes.
                    as.integer(rs$riskset),
                    #
                    as.integer(nn),
                    #
                    as.double(score),       ## 'score'
                    as.double(hazard),
                    resid = double(nn)
                    #DUP = FALSE
                    )$resid

  if (!fit$fail)
    var <- matrix(fit$variance, ncov, ncov)
  else
    var <- NULL

  bootstrap <- NULL
  boot.sd <- NULL
  if (boot & (fit$fail == 0)){
    if (!is.numeric(boot)){
      cat("boot must be numeric (number of bootstrap replicates)")
    }else{
      init <- fit$beta
      iter <- control$maxiter
      fit.boot <- .Fortran("bootcox",
                           as.integer(1), ## means 'coxreg'
                           as.integer(boot),
                           boot.sample = double(boot * ncov),
                           boot.sd = double(boot * ncov),
                           as.integer(method == "efron"),
                           iter = as.integer(iter),
                           as.double(control$eps),
                           as.integer(printlevel),
                                        #
                           as.integer(sum(rs$n.events)), 
                           as.integer(sum(rs$antrs)),  
                           as.integer(length(rs$antrs)),
                                        #
                           as.integer(rs$antrs),
                           as.integer(rs$n.events),
                           as.integer(rs$size),
                                        #
                           as.integer(length(rs$riskset)), 
                           as.integer(rs$eventset),
                           as.integer(rs$riskset),
                                        #
                           as.integer(nn),
                           as.integer(ncov),
                           as.double(scale(X, center = TRUE, scale = FALSE)),
                           as.double(offset),
                                        #
                           as.double(init),     
                           as.double(fit$beta), ## Estimated beta
                                        #
                           loglik = double(2), 
                           dloglik = double(ncov),
                           variance = double(ncov * ncov),
                           sctest = double(1),
                                        #
                           double(nn),     
                           double(ncov),   
                           double(ncov * ncov),
                                        #
                           conver = integer(1),
                           fail = integer(1)
                           ##DUP = FALSE,
                           )
      bootstrap <- matrix(fit.boot$boot.sample, ncol = ncov, byrow = TRUE)
      boot.sd <- matrix(fit.boot$boot.sd, ncol = ncov, byrow = TRUE)
    }      
   }   

  list(coefficients = fit$beta,
       var = var,
       loglik = fit$loglik,
       score = fit$sctest,
       linear.predictors = X %*% fit$beta,
       residuals = resid,
       hazard = hazard,
       means = apply(X, 2, mean),
       bootstrap = bootstrap,
       boot.sd = boot.sd,
       conver = fit$conver,
       fail = fit$fail,
       iter = fit$iter
       )
       
  }
