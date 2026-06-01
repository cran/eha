#' Parametric AFT regression
#' 
#' This function is called by \code{\link{aftreg}}, but it can also be directly
#' called by a user.
#' 
#' See \code{\link{aftreg}} for more detail.
#' 
#' @param X The design (covariate) matrix.
#' @param Y A survival object, the response.
#' @param param Which parametrization?
#' @param strata A stratum variable.
#' @param offset Offset.
#' @param init Initial regression parameter values.
#' @param id See corresponding argument to \code{\link{aftreg}}.
#' @param control Controls convergence and output.
#' @return 
#' \item{coefficients}{Estimated regression coefficients plus estimated
#' scale and shape coefficients, sorted by strata, if present.}
#' \item{df}{Degrees of freedom; No. of regression parameters.}
#' \item{var}{Variance-covariance matrix} 
#' \item{loglik}{Vector of length 2. The
#' first component is the maximized loglihood with only scale and shape in the
#' model, the second the final maximum.} 
#' \item{conver}{TRUE if convergence} 
#' \item{fail}{TRUE if failure} 
#' \item{iter}{Number of Newton-Raphson iterates.} 
#' \item{n.strata}{The number of strata in the data.}
#' @author Göran Broström
#' @seealso \code{\link{aftreg}}
#' @keywords survival regression
#' @export
aftGompreg.fit <- function(X, Y, param,
                           strata, offset,
                           init, id,
                           control){

    ## New in Version 1.2-9; wrong before!
    ## Note that we MUST keep individuals together here;
    ## stratum comes second, beacuse we can then let individuals
    ## change stratum over time!
    ## cat("length(id) = ", length(id), "dim(Y) = ", dim(Y), "\n")
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
    dist <- "gompertz"
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


    if (dist == "gompertz"){ # Special treatment
        pfixed <- FALSE # No fix for Gompertz (yet).
        fit <- aftGomp0(printlevel, ns, nn, id,
                      strata, Y, X, offset, means, param)
    }            

    ## Add back means:
    if (ncov){
        rew <- addMeans(means = means,
                        par = fit$beta,
                        var = fit$var,
                        ns = ns,
                        pfixed = FALSE, # Never TRUE for Gompertz!
                        coef.names = colnames(X))
        fit$beta <- rew$par
        fit$var <- rew$var
    }else{
        ## Just add names and scale, shape
        if (ns > 1){
            coef.names <- c("log(scale):1")
            coef.names <- c(coef.names, "log(shape):1")
            for (i in 2:ns){
                coef.names <- c(coef.names,
                                paste("log(scale)", as.character(i), sep =":"))
                coef.names <- c(coef.names, 
                                paste("log(shape)", as.character(i), sep =":"))
            }
            
        }else{
            coef.names <- "log(scale)"
            coef.names <- c(coef.names, "log(shape)")
        }
        names(fit$beta) <- coef.names
        colnames(fit$var) <- rownames(fit$var) <- coef.names
    }

    ## Reparametrisation (if asked for):

    if ((param == "lifeExp") && ncov){
        new.par <- translifeExp(fit$beta, fit$var, ns, pfixed)
        coefficients <- new.par$coefficients
        var <- new.par$var
    }else{
        coefficients <- fit$beta
        var <- fit$var
    }
    ##cat("done!\n")
    scale <- numeric(ns)
    shape <- numeric(ns)
    j <- 0
    for (i in 1:ns){
        j <- j + 1
        scale[i] <- exp(coefficients[ncov + j])
        j <- j + 1
        shape[i] <- exp(coefficients[ncov + j])
    }
    ##coefficients <- fit$beta

    list(coefficients = coefficients,
         df = fit$ncov,
         var = var,
         loglik = fit$loglik,
         convergence = !fit$fail,
         fail = fit$fail,
         n.strata = ns,
         scale = scale,
         shape = shape
         )
}
