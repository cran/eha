weibreg.fit <- function(X, Y, strata, offset, init, shape, control){

  nn <- NROW(X)
  ncov <- NCOL(X)

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

  if (shape <= 0){

      bdim <- ncov + 2 * ns
      if (ns > 0){
          ord <- order(strata)
          X <- X[ord, ]
          Y <- Y[ord, ]
          offset <- offset[ord]
          nstra <- c(0, cumsum(table(strata)))
      }
      X <- scale(X, center = TRUE, scale = FALSE)
      fit <- .C("weibsup",
                iter = as.integer(iter), #maxit on ip, actual on op.
                as.double(control$eps),
                as.integer(printlevel),
                                        #
                as.integer(ns), # No. of strata
                as.integer(nstra),
                as.integer(nn),
                as.integer(ncov),
                as.integer(bdim),
                                        #
                as.double(Y[, 1]),  ## 'enter'
                as.double(Y[, 2]),  ## 'exit'
                as.integer(Y[, 3]), ## 'event'
                                        #
                as.double(t(X)), ## NOTE transpose!
                as.double(offset),
                                        #
                as.double(init),     # 'start.beta'
                beta = double(bdim), # results -->
                lambda = double(ns),
                lambda.sd = double(ns),
                shape = double(ns),  ## "p"
                shape.sd = double(ns),
                                        #
                loglik = double(2), # [1] == start, [2] == maximized
                dloglik = double(bdim),
                variance = double(bdim * bdim),
                sctest = double(1),
                                        #
                conver = integer(1),
                fail = integer(1),
                #DUP = FALSE,
                PACKAGE = "eha")

      if (fit$fail) return(fail = fit$fail,
                           n.strata = ns,
                           value = fit$beta[fit$fail])
      
      for (i in 1:ns) ## Really a HACK !!!!!!!!!!!!!!!
          fit$beta[ncov + 2 * i - 1] <- -fit$beta[ncov + 2 * i - 1]
      fit$shape.fixed <- FALSE
  }else{  ## Exponential regression:
      if (ns >= 2) warning("'strata' is not meaningful for exponential regression.\n Include stratum variable as a factor in the model instead.")
      bdim <- ncov + 1
      X <- scale(X, center = TRUE, scale = FALSE)

      fit <- .C("expsup",
                iter = as.integer(iter), #maxit on ip, actual on op.
                as.double(control$eps),
                as.integer(printlevel),
                                        #
                as.integer(nn),
                as.integer(ncov),
                as.integer(bdim),
                #
                as.double(Y[, 1]),  ## 'enter'
                as.double(Y[, 2]),  ## 'exit'
                as.integer(Y[, 3]), ## 'event'
                #
                as.double(t(X)),
                as.double(offset),
                as.double(shape), ## "p"
                                        #
                as.double(init),     # 'start.beta'
                beta = double(bdim), # results -->
                lambda = double(1),
                lambda.sd = double(1),
                                        #
                loglik = double(2), # [1] == start, [2] == maximized
                dloglik = double(bdim),
                variance = double(bdim * bdim),
                sctest = double(1),
                                        #
                conver = integer(1),
                fail = integer(1),
                                        #DUP = FALSE,
                PACKAGE = "eha")
      if (fit$fail) return(fail = fit$fail,
                           n.strata = ns,
                           value = fit$beta[fit$fail])
      fit$shape.fixed <- TRUE
      fit$shape <- shape
      fit$shape.sd <- NULL  ## Not necessary!?!?
      fit$beta[bdim] <- -fit$beta[bdim] ## To get "1 / lambda"!
      ## Note; this is really a "hack"!!!!!!!!!!!!!!!
}

  lp <- offset + X %*% fit$beta[1:ncov]
  score <- exp(lp)
  ##cat("done!\n")
      
  if (!fit$fail)
    var <- matrix(fit$variance, bdim, bdim)
  else
    var <- NULL

  
  list(coefficients = fit$beta,
       var = var,
       loglik = fit$loglik,
       score = fit$sctest,
       linear.predictors = lp,
       means = apply(X, 2, mean),
       conver = fit$conver,
       fail = fit$fail,
       iter = fit$iter,
       n.strata = ns
       )
       
}
