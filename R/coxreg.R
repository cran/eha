coxreg <- function (formula = formula(data),
                    data = parent.frame(), 
                    na.action = getOption("na.action"),
                    init = NULL,
                    method = c("efron", "breslow", "mppl", "ml"),
                    control = list(eps = 1e-8, maxiter = 25, trace = FALSE),
                    singular.ok = TRUE,
                    model = FALSE, 
                    x = FALSE,
                    y = TRUE,
                    boot = FALSE,
                    geometric = NULL,
                    rs = NULL,
                    frailty = NULL,
                    max.survs = NULL) 
{
    if (!is.null(frailty))
      stop("Frailty not implemented yet.")
    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
 
    special <- "strata"
    Terms <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    if (is.null(max.survs)) max.survs <- NROW(Y)
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
        temp <- untangle.specials(Terms, "strata", 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars) == 1) 
            strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
    if (length(dropx)) 
        newTerms <- Terms[-dropx]
    else newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 
        1)
    X <- X[, -1, drop = FALSE]

    #########################################

    if (length(dropx)){
      covars <- names(m)[-c(1, (dropx + 1))]
    }else{
      covars <- names(m)[-1]
    }

    isF <- logical(length(covars))
    if (length(covars)){
        for (i in 1:length(covars)){
            if (length(dropx)){
                isF[i] <- ( is.factor(m[, -(dropx + 1)][, (i + 1)]) ||
                           is.logical(m[, -(dropx + 1)][, (i + 1)]) )
            }else{
                isF[i] <- ( is.factor(m[, (i + 1)]) ||
                           is.logical(m[, (i + 1)]) )
            }      
        }
    }
    
    if (any(isF)){
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

    ##########################################
    type <- attr(Y, "type")
    if (type != "right" && type != "counting") 
      stop(paste("Cox model doesn't support \"", type, "\" survival data", 
                 sep = ""))
    
    if (NCOL(Y) == 2){
        Y <- cbind(numeric(NROW(Y)), Y)
    }
    
    n.events <- sum(Y[, 3] != 0)
    if (n.events == 0) stop("No events; no sense in continuing!")

    if ((!is.null(init)) && (length(init) != NCOL(X)))
      stop("Wrong length of 'init'")
    
    
    if (is.list(control)){
        if (is.null(control$eps)) control$eps <- 1e-8
        if (is.null(control$maxiter)) control$maxiter <- 10
        if (is.null(control$trace)) control$trace <- FALSE
    }else{
        stop("control must be a list")
    }
    
    fit <- coxreg.fit(X,
                      Y,
                      rs,
                      strats,
                      offset,
                      init,
                      max.survs,
                      method,
                      boot,
                      calc.hazards = TRUE,
                      calc.martres = TRUE,
                      control,
                      verbose = TRUE)

##    if (!length(fit$coefficients)){
##        class(fit) <- c("coxreg", "coxph")
##        return(fit)
##    }
    ##if (is.null(fit)) return(NULL) ## Removed 19 Feb 2007
    ##if (!fit$fail) fit$fail <- NULL
    ##else
    ##    fit$fail <- TRUE

    fit$convergence <- as.logical(fit$conver)
    fit$conver <- NULL ## Ugly!
    fit$f.convergence <- as.logical(fit$f.conver)
    fit$f.conver <- NULL
###########################################################################    
## Crap dealt with ......
    
    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "coxreg"
    }
    else if (fit$fail){
        if (length(fit$coef) && any(is.na(fit$coef))) {
            vars <- (1:length(fit$coef))[is.na(fit$coef)]
            msg <- paste("X matrix deemed to be singular; variable", 
                paste(vars, collapse = " "))
            if (singular.ok) 
                warning(msg)
            else stop(msg)
        }
        fit$n <- nrow(Y)
        class(fit) <- fit$method
        fit$terms <- Terms
        fit$assign <- assign
        if (FALSE){ ## Out-commented
        if (length(fit$coef) && is.null(fit$wald.test)) {
            nabeta <- !is.na(fit$coef)
            if (is.null(init)) 
                temp <- fit$coef[nabeta]
            else temp <- (fit$coef - init)[nabeta]
            fit$wald.test <- survival:::coxph.wtest(fit$var[nabeta, nabeta], 
                temp, control$toler.chol)$test
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

    fit$isF <- isF
    fit$covars <- covars
    s.wght <- (Y[, 2] - Y[, 1])## * weights
    fit$ttr <- sum(s.wght)
    fit$w.means <- list()
    if (length(fit$covars)){
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
                fit$w.means[[i]] <- sum(s.wght * m[, col.m]) / fit$ttr
            }
        }
    }
    ##########################################
    fit$levels <- levels
    fit$formula <- formula(Terms)
    fit$terms <- Terms
    fit$call <- call
    fit$events <- n.events
    if (length(fit$coefficients)){
        names(fit$coefficients) <- colnames(X)
        fit$means <- apply(X, 2, mean)
    }
    fit$method <- method
    class(fit) <- c("coxreg", "coxph")
    fit
}
