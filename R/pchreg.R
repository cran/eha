#' Piecewise Constant Proportional Hazards Regression
#' 
#' Proportional hazards model with piecewise constant baseline hazard(s).  Allows for
#' stratification and left truncated and right censored data.
#' 
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right.  The response must be a survival
#' object as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param cuts Specifies the points in time
#' where the hazard function jumps. If omitted, an exponential model is fitted.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param control a list with components \code{eps} (convergence criterion),
#' \code{maxiter} (maximum number of iterations), and \code{silent} (logical,
#' controlling amount of output).  You can change any component without mention
#' the other(s).
#' @param singular.ok Not used.
#' @param model Not used.
#' @param x Return the design matrix in the model object?
#' @param y Return the response in the model object?
#' @return A list of class \code{"pchreg"} with components
#' \item{coefficients}{Fitted parameter estimates.} 
#' \item{cuts}{Cut points (\code{NULL} if no cut points).} 
#' \item{hazards}{The estimated constant levels.} 
#' \item{var}{Covariance matrix of the estimates.}
#' \item{loglik}{Vector of length two; first component is the value at the
#' initial parameter values, the second component is the maximized value.}
#' \item{score}{The score test statistic (at the initial value).}
#' \item{linear.predictors}{The estimated linear predictors.}
#' \item{means}{Means of the columns of the design matrix, except those columns
#' corresponding to a factor level. Otherwise all
#' zero.} 
#' \item{w.means}{Weighted (against exposure time) means of covariates;
#' weighted relative frequencies of levels of factors.} 
#' \item{n}{Number of
#' spells in indata (possibly after removal of cases with NA's).}
#' \item{n.events}{Number of events in data.} 
#' \item{terms}{Used by extractor functions.} 
#' \item{assign}{Used by extractor functions.} %
#' \item{wald.test}{The Wald test statistic (at the initial value).}
#' \item{y}{The Surv vector.} 
#' \item{isF}{Logical vector indicating the
#' covariates that are factors.} 
#' \item{covars}{The covariates.}
#' \item{ttr}{Total Time at Risk.} 
#' \item{levels}{List of levels of factors.}
#' \item{formula}{The calling formula.} 
#' \item{call}{The call.}
#' \item{method}{The method.} 
#' \item{convergence}{Did the optimization
#' converge?} 
#' \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
#' @author Göran Broström
#' @seealso \code{\link{phreg}}, \code{\link{coxreg}}, \code{link{aftreg}}.
#' @keywords survival regression
#' @examples
#' 
#' \dontrun{
#' dat <- age.window(oldmort, c(60, 80))
#' fit <- pchreg(Surv(enter, exit, event) ~ ses.50 + sex, 
#' data = dat, cuts = seq(60, 80, by = 4))
#' summary(fit)
#' 
#' fit.cr <- coxreg(Surv(enter, exit, event) ~ ses.50 + sex, data = dat)
#' check.dist(fit.cr, fit, main = "Cumulative hazards")
#' }
#' 
#' @export pchreg
pchreg <- function (formula = formula(data),
                    data = parent.frame(),
                    na.action = getOption("na.action"),
                    cuts = NULL,
                    init,
                    control = list(eps = 1e-8, maxiter = 20, trace = FALSE),
                    singular.ok = TRUE,
                    model = FALSE,
                    x = FALSE,
                    y = TRUE)
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  temp <- c("", "formula", "data", "na.action")
  m <- m[match(temp, names(m), nomatch = 0)] # m is a call
  
  special <- "strata"
  Terms <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data) # Terms is a 'terms' 'formula'
  
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  ##return(m)
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  
  weights <- model.extract(m, "weights")
  offset <- attr(Terms, "offset")
  tt <- length(offset)
  offset <- if (tt == 0)
    rep(0, nrow(Y))
  else if (tt == 1)
    m[[offset]]
  else {
    ff <- m[[offset[1]]]
    for (i in 2:tt) ff <- ff + m[[offset[i]]]
    ff
  }
  attr(Terms, "intercept") <- 1
  strats <- attr(Terms, "specials")$strata
  dropx <- NULL
  
  if (length(strats)) {
    ##if (dist == "pch") # Changed 2.4-0
    ##  stop("No strata allowed in the pch model (yet)") 
    temp <- survival::untangle.specials(Terms, "strata", 1)
    dropx <- c(dropx, temp$terms)
    if (length(temp$vars) == 1)
      strata.keep <- m[[temp$vars]]
    else strata.keep <- survival::strata(m[, temp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
  }
  
  if (length(dropx))
    newTerms <- Terms[-dropx]
  else newTerms <- Terms
  X <- model.matrix(newTerms, m)
  ##return(X)
  assign <- lapply(survival::attrassign(X, newTerms), function(x) x - 1)
  
  X <- X[, -1, drop = FALSE]
  intercept <- FALSE
  ncov <- NCOL(X)
  nullModel <- ncov == 0
  #########################################
  
  if (ncov){
    if (length(dropx)){
      covars <- names(m)[-c(1, (dropx + 1))]
    }else{
      covars <- names(m)[-1]
    }
    
    isF <- logical(length(covars))
    for (i in 1:length(covars)){
      if (length(dropx)){
        if (is.logical(m[, -(dropx + 1)][, (i + 1)])){
          m[, -(dropx + 1)][, (i + 1)] <-
            as.factor(m[, -(dropx + 1)][, (i + 1)])
        }
        
        isF[i] <- is.factor(m[, -(dropx + 1)][, (i + 1)])## ||
        ##is.logical(m[, -(dropx + 1)][, (i + 1)]) )
      }else{
        if (is.logical(m[, (i + 1)])){
          m[, (i + 1)] <- as.factor(m[, (i + 1)])
        }
        
        isF[i] <- is.factor(m[, (i + 1)])## ||
        ##is.logical(m[, (i + 1)]) )
      }
    }
    
    ant.fak <- sum(isF)
  }else{ #!ncov
    isF <- logical(0)
    ant.fak <- 0
  }
  ##cat("ant.fak = ", ant.fak, "\n")
  if (ant.fak){
    levels <- list()
    index <- 0
    for ( i in 1:length(covars) ){
      if (isF[i]){
        index <- index + 1
        if (length(dropx)){
          levels[[i]] <- levels(m[, -(dropx + 1)][, (i + 1)])
        }else{
          levels[[i]] <- levels(m[, (i + 1)])
        }
      }else{
        levels[[i]] <- NULL
      }
    }
  }else{
    levels <- NULL
  }
  
  isI <- logical(NCOL(X))
  if (ant.fak){
    indx <- 0
    for (i in seq_len(length(covars))){
      indx <- indx + 1
      if (isF[i]){
        isI[indx] <- TRUE
        isI[indx] <- TRUE
        if (length(levels[[i]]) >= 3){
          for (j in 3:length(levels[[i]])){
            indx <- indx + 1
            isI[indx] <- TRUE
          }
        }
      }
    }
  }

  ##########################################
  type <- attr(Y, "type")
  if (type != "right" && type != "counting")
    stop(paste("This model doesn't support \"", type, "\" survival data",
               sep = ""))
  
  if (NCOL(Y) == 2){
    Y <- cbind(numeric(NROW(Y)), Y)
  }
  
  n.events <- sum(Y[, 3] != 0)
  if (n.events == 0) stop("No events; no sense in continuing!")
  if (missing(init)){ # Is this wise? No!
    if (ncov){
      ##init <- coxreg(formula, data = data)$coefficients
      init <- rep(0, ncov)
    }else{
      init <- numeric(0)
    }
  }
  if (is.list(control)){
    if (is.null(control$eps)) control$eps <- 1e-8
    if (is.null(control$maxiter)) control$maxiter <- 10
    if (is.null(control$trace)) control$trace <- FALSE
  }else{
    stop("control must be a list")
  }
  
  if (missing(cuts)){
    ##stop("'dist = pch' needs 'cuts' to be set")
    cuts <- numeric(0) # Exponential distribution(s)
  }
  fit <- pchreg.fit(X,
                    Y,
                    cuts,
                    offset,
                    strats,
                    init,
                    control)
  if (fit$fail){
    warning(paste("Failed with error code ", fit$fail))
    return(1)
  }
  
  if (ncov){
    ##cat("ncov = ", ncov, "\n")
    fit$linear.predictors <- offset + X %*%
      fit$coefficients[1:ncov]
  }else{
    fit$linear.predictors <- numeric(0)
    ##fit$means <- numeric(0)
  }
  ##score <- exp(lp)
  
  ##cat("fit$means == ", fit$means, "\n")
  if (!fit$fail){
    fit$fail <- NULL
  }else{
    out <- paste("Singular hessian; suspicious variable No. ",
                 as.character(fit$fail), ":\n",
                 names(coefficients)[fit$fail], " = ",
                 as.character(fit$value),
                 "\nTry running with fixed shape", sep = "")
    stop(out)
  }
  
  
  fit$convergence <- as.logical(fit$conver)
  fit$conver <- NULL ## Ugly!
  ###########################################################################
  ## Crap dealt with ......
  
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "mlreg"
  }
  else if (is.null(fit$fail)){
    if (!is.null(fit$coef) && any(is.na(fit$coef))) {
      vars <- (1:length(fit$coef))[is.na(fit$coef)]
      msg <- paste("X matrix deemed to be singular; variable",
                   paste(vars, collapse = " "))
      if (singular.ok)
        warning(msg)
      else stop(msg)
    }
    fit$n <- nrow(Y)
    fit$terms <- Terms
    fit$assign <- assign
    if (FALSE){ ## Out-commented...(why?)
      if (length(fit$coef) && is.null(fit$wald.test)) {
        nabeta <- !is.na(fit$coef)
        if (is.null(init))
          temp <- fit$coef[nabeta]
        else temp <- (fit$coef - init)[nabeta]
        ##fit$wald.test <-
        ##survival:::coxph.wtest(fit$var[nabeta, nabeta],
        ##                     temp, control$toler.chol)$test
      }
    }
    na.action <- attr(m, "na.action")
    if (length(na.action))
      fit$na.action <- na.action
    if (model)
      fit$model <- m
    if (x) {
      fit$x <- X
      if (length(strats))
        fit$strata <- strata.keep
    }
    if (y)
      fit$y <- Y
  }
  ##if (!is.null(weights) && any(weights != 1))
  ##    fit$weights <- weights
  
  ##########################################
  s.wght <- (Y[, 2] - Y[, 1])## * weights
  fit$ttr <- sum(s.wght)
  if (ncov){
    fit$isF <- isF
    fit$covars <- covars
    fit$w.means <- list()
    for (i in 1:length(fit$covars)){
      nam <- fit$covars[i]
      col.m <- which(nam == names(m))
      if (isF[i]){
        n.lev <- length(levels[[i]])
        fit$w.means[[i]] <- numeric(n.lev)
        for (j in 1:n.lev){
          who <- m[, col.m] == levels[[i]][j]
          fit$w.means[[i]][j] <-
            sum( s.wght[who] ) / fit$ttr ## * 100, if in per cent
        }
      }else{
        ##fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
        fit$w.means[[i]] <- weighted.mean(m[, col.m], s.wght)
      }
    }
  }
  
  ##########################################
  fit$ttr <- sum(s.wght)
  ##names(fit$coefficients) <- coef.names
  fit$levels <- levels
  fit$formula <- formula(Terms)
  fit$call <- call
  fit$dist <- "pch"
  fit$n.events <- n.events
  fit$nullModel <- nullModel
  class(fit) <- c("pchreg", "tpchreg")
  if (length(strats))
    fit$strata <- names(strats)
  if (length(strats)){
    fit$strata <- levels(as.factor(strata.keep))
  }
  fit
}
#' @export
extractAIC.pchreg <- function(fit, scale, k = 2, ...){
  edf <- sum(fit$df) + length(fit$hazards)
  loglik <- fit$loglik[length(fit$loglik)]
  c(edf, -2 * loglik + k * edf)
}

#' @export
nobs.pchreg <- function(object, ...){
  object$n
}

