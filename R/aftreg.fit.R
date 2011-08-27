aftreg.fit <- function(X, Y, dist,
                       strata, offset,
                       init, shape, id,
                       control, center = NULL){

    ## New in Version 1.2-9; wrong before!
    ## Note that we MUST keep individuals together here;
    ## stratum comes second, beacuse we can then let individuals
    ## change stratum over time!
    ord <- order(id, Y[, 1])
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    if (length(colnames(X)) != NCOL(X)){
        colnames(X) <- paste("X", 1:NCOL(X), sep = ".")
    }
    if (NROW(X) != NROW(Y)) stop("Wrong No. of rows in X")
    X <- X[ord, , drop = FALSE]
    Y <- Y[ord, , drop = FALSE]
    id <- id[ord]
    strata <- strata[ord]
    offset <- offset[ord]

##### Check 'data integrity' ########
    stop.flag <- FALSE
    overlap <- NULL
    n.ind <- length(unique(id))
    if (length(id) > NROW(Y)) stop("'id' argument too long")
    if (any(Y[, 2] <= Y[, 1])) stop("Zero or negative length intervals")
    else if (n.ind < NROW(Y)){
        for (i in 2:NROW(Y)){
            if (id[i-1] == id[i]){
                if (Y[i-1, 2] > Y[i, 1]){
                    cat(paste("Overlapping intervals for id ", id[i], "\n"))
                    overlap <- c(overlap, id[i])
                    stop.flag <- TRUE
                }
            }
        }
    }
    if (stop.flag) {
        cat("Error(s) in data encountered\n")
        fit <- list()
        fit$overlap <- overlap
        return(fit)
    }
#####################################
    if (dist == "weibull"){
        dis <- 0
    }else if(dist == "loglogistic"){
        dis <- 1
    }else if (dist == "lognormal"){
        dis <- 2
    }else if (dist == "ev"){
        dis <- 3
    }else if (dist == "gompertz"){ # An EV with shape == 1: NOT REALLY!!
        dis <- 4
        ##dis <- 3
        ##shape <- 1
        ##stop("The Gompertz is not available yet as AFT; try PH?")
    }else{
        stop(paste(dist, "is not an implemented distribution"))
    }

    nn <- NROW(X)
    ncov <- NCOL(X) ## No intercept! - 1

    ## No intercepts in 'aftreg!!! (1.2-17) intercept <- (dis == 4) # gompertz
    ## Again: No intercept in X:
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

    if (all(shape <= 0) || (dist == "gompertz")){
        ## Then shape is estimated in all strata
        if (dist == "gompertz"){
            fit <- aftp0g(printlevel, ns, nn, id,
                         strata, Y, X, offset, means)
        }else{
            fit <- aftp0(printlevel, ns, nn, id,
                         strata, Y, X, offset, dis, means)
        }
    }else{
        ## Then shape is fixed in all strata:
        ## Note: We must allow stratification even here (091006)!!
        fit <- aftp1(printlevel, ns, nn, id,
                     strata, Y, X, offset, shape, dis, means)
    }            
    

    ##cat("done!\n")

    coefficients <- fit$beta
    names(coefficients) <- fit$coef.names

    list(coefficients = coefficients,
         df = fit$ncov,
         var = fit$var,
         loglik = fit$loglik,
         ##score = fit$sctest,
         convergence = (fit$conver == 0),
         fail = fit$fail,
         ##iter = fit$iter,
         n.strata = ns#,
         ##shape = fit$shape
         )

}
