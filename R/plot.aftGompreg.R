#' Plots output from an AFT regression, aftGompreg.
#' 
#' Just a simple plot of the hazard (cumulative hazard, density, survival)
#' functions for each stratum.
#' 
#' The plot is drawn at the mean values of the covariates, by default.
#' 
#' @param x A \code{aftGompreg} object
#' @param fn Which functions shoud be plotted! Default is the hazard function(s)
#' @param main Header for the plot
#' @param xlim x limits
#' @param ylim y limits
#' @param xlab x label
#' @param ylab y label
#' @param col Colors?
#' @param lty Line types?
#' @param printLegend Logical, or character ("topleft", "bottomleft",
#' "topright" or "bottomright"); if \code{TRUE} or character, a legend is added
#' to the plot if the number of strata is two or more.
#' @param new.data At which covariate values?
#' @param \dots Extra parameters passed to 'plot'
#' @return No return value.
#' @author Göran Broström
#' @seealso \code{\link{aftreg}}
#' @keywords dplot survival
#' @examples
#' 
#' y <- rllogis(40, shape = 1, scale = 1)
#' x <- rep(c(1,1,2,2), 10)
#' fit <- aftreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
#' plot(fit)
#' 
#' @export
plot.aftGompreg <- function(x,
                        fn = c("haz", "cum", "den", "sur"),
                        main = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        xlab = "Duration",
                        ylab = "",
                        col,
                        lty,
                        printLegend = TRUE,
                        new.data = x$means,
                         ...){
    if (!inherits(x, "aftGompreg")) stop("Works only with 'aftGompreg' objects.")
    ##if (x$pfixed) stop("True exponential hazards are not plotted")
    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if (missing(col)) col <- rep(1, x$n.strata) ## New 2013-12-05
    if (missing(lty)) lty <- 1:x$n.strata # No. of strata

    if (length(col) < x$n.strata) col <- rep(col, x$n.strata)
    if (length(lty) < x$n.strata) lty <- rep(lty, x$n.strata)

 #   fn <- fn[1]

    if (length(fn) >= 3){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }else if (length(fn) == 2){
        oldpar <- par(mfrow = c(2, 1))
        on.exit(par(oldpar))
    }
    ##ncov <- length(x$means) # Doesn't work with some factor covariates!!
    ncov <- x$df
    ns <- x$n.strata
    param.scale <- if (x$param=="lifeAcc") -1 else 1
### Check:
    p <- exp(x$coefficients[ncov + (1:ns) * 2])
    lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1] +
                  param.scale*sum((new.data - x$means) * x$coefficients[1:ncov]))
###

    if (ncov){
        score <- exp(sum((new.data - x$means) * x$coefficients[1:ncov]))
    }else{
        score <- 1
    }

    ##if (ncov){ # THIS IS for aftplot!!
    ##    uppe <- exp(-sum((new.data[1:ncov] * x$coefficients[1:ncov]) / p)
    ##    lambda <- lambda * uppe
    ##}
    if (is.null(xlim))
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))

    npts <- 4999
    xx <- seq(xlim[1], xlim[2], length = npts)
    ##if (xx[1] <= 0) xx[1] <- 0.001

    ## hazard
 if (x$dist == "gompertz"){
        dist = "Gompertz"
        haza <- hgompertz
        Haza <- Hgompertz
        Surviv <- pgompertz
        Dens <- dgompertz
        ## canonical parameterisation
        if (x$param == "canonical") p <- p/lambda
 }
    ##cat("shape = ", p, ", scale = ", lambda, "\n")

    if ("haz" %in% fn){
        haz <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
            haz[i, ] <- haza(xx, scale = lambda[i], shape = p[i],
                             param = "canonical")
        }

        if (is.null(ylim)) {
            ylim <- c(0, max(haz))
            if (min(p) < 1) ylim[2] <- min(ylim[2], max(haz[, -1]))
        }
        ##if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)) main <- paste(dist, "hazard function")
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main,
             col = col, lty = lty, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = lty[i], col = col[i], ...)
            }
        }
        abline(h = 0)
        abline(v = 0)

    }
    ## Cumulative hazard
    if ("cum" %in% fn){

        Haz <- matrix(ncol = npts, nrow = ns)

    ##
        for (i in 1:ns){
            Haz[i, ] <- Haza(xx, scale = lambda[i], shape = p[i],
                             param = "canonical")
        }
        if (is.null(ylim)) ylim <- c(0, max(Haz))
        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazards"
        if (is.null(main))
        main <- paste(dist, "cumulative hazard function")
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, col = col, lty = lty, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = lty[i], col = col[i], ...)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## density
    if ("den" %in% fn){

        den <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
           den[i, ] <- Dens(xx, scale = lambda[i], shape = p[i], 
                            param = "canonical")
        }

        if (is.null(ylim)){
            ylim <- c(0, max(den))
            if (min(p) < 1) ylim[2] <- min(max(den[, -1]))
        }
        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Density"
        if (is.null(main))
            main <- paste(dist, "density function")
        plot(xx, den[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, col = col, lty = lty, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, den[i, ], type = "l", lty = lty[i], col = col[i], ...)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    ## Survivor function
    if ("sur" %in% fn){


        sur <- matrix(ncol = npts, nrow = ns)
        for (i in 1:ns){
                sur[i, ] <- Surviv(xx, scale = lambda[i],
                                   shape = p[i], param = "canonical",
                                   lower.tail = FALSE)
            
        }

        if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        if (is.null(main))
            main <- paste(dist, "survivor function")
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, col = col, lty = lty, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = lty[i], col = col[i], ...)
            }
        }
        abline(h = 0)
        abline(v = 0)
    }
    
    if (is.character(printLegend)){
        if (!(printLegend %in% c("topleft", "bottomleft",
                                 "topright", "bottomright",
                                 "bottom", "left", "top",
                                 "right", "center"))){
            printLegend <- FALSE
            warning("Illegal value of 'printLegend'")
        }
    }
    if (is.logical(printLegend)){
        if ((ns > 1) && printLegend){
            legend(x = "bottomright",  legend = x$strata, lty = lty,
                   inset = 0.001,
                   col = col)
        }
    }else{
        if ((ns > 1) && is.character(printLegend)){
            legend(x = printLegend,  legend = x$strata,
                   lty = lty, inset = 0.001, col = col)
        }
    }
    
    ##par(oldpar)
}
