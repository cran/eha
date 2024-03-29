#' Plots output from a phreg regression
#' 
#' Plot(s) of the hazard, density, cumulative hazards, and/or the survivor
#' function(s) for each stratum.
#' 
#' 
#' @param x A \code{phreg} object
#' @param fn Which function should be plotted? Default is the hazard function(s).  
#' @param main Header for the plot
#' @param xlim x limits
#' @param ylim y limits
#' @param xlab x label
#' @param ylab y label
#' @param col Color(s) for the curves. Defaults to black.
#' @param lty Line type for the curve(s). Defaults to 1:(No. of strata).
#' @param printLegend Logical, or character ("topleft", "bottomleft",
#' "topright" or "bottomright"); if \code{TRUE} or character, a legend is added
#' to the plot if the number of strata is two or more.
#' @param score Multiplication factor for the hazard function.
#' @param fig logical, should the graph be drawn? If FALSE, data is returned.
#' @param \dots Extra parameters passed to 'plot' and 'lines'.
#' @return No return value if fig = TRUE, otherwise the cumulative 
#' hazards function (coordinates), given \code{fn = "cum"}.
#' @note Reference hazard is given by the fit; 
#' zero for all covariates, and  the reference category for factors.
#' @author Göran Broström
#' @seealso \code{\link{phreg}}
#' @keywords dplot survival
#' @examples
#' 
#' y <- rllogis(40, shape = 1, scale = 1)
#' x <- rep(c(1,1,2,2), 10)
#' fit <- phreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
#' plot(fit)
#' 
#' @export
plot.phreg <- function(x,
                       fn = c("haz", "cum", "den", "sur"),
                       main = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       xlab = "Duration",
                       ylab = "",
                       col,   ## New 6 Feb 2013
                       lty,   ## New 6 Feb 2013
                       printLegend = TRUE,
                       ##legend = printLegend,
                       score = 1,
                       fig = TRUE,
                       ...){

    if (!inherits(x, "phreg")) stop("Works only with 'phreg' objects.")

    if (missing(col)) col <- rep(1, x$n.strata) ## New 2013-12-05
    if (missing(lty)) lty <- 1:x$n.strata # No. of strata

    if (length(col) < x$n.strata) col <- rep(col, x$n.strata)
    if (length(lty) < x$n.strata) lty <- rep(lty, x$n.strata)

    fn <- fn[1]
    if (fn == "surv") fn <- "sur"
    if (fn == "cumhaz") fn <- "cum"
    if (!(fn %in% c("haz", "cum", "den", "sur")))
        stop(paste(fn, "is an illegal value of 'fn'"))


    ##if (length(fn) >= 3){
    ##    oldpar <- par(mfrow = c(2, 2))
    ##    on.exit(par(oldpar))
    ##}else if (length(fn) == 2){
    ##    oldpar <- par(mfrow = c(2, 1))
    ##    on.exit(par(oldpar))
    ##}
    ##ncov <- length(x$w.means) # Doesn't work with some factor covariates!!
    ncov <- x$df # Better!
    ns <- x$n.strata
    if (!is.null(x$pfixed)){
        if (x$pfixed){
            shape <- rep(x$shape, ns)
            scale <- exp(x$coefficients[ncov + (1:ns)])
        }else{
            shape <- exp(x$coefficients[ncov + (1:ns) * 2])
            scale <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
        }
    }else if (x$dist != "pch"){
        shape <- exp(x$coefficients[ncov + (1:ns) * 2])
        scale <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    }


    ##if (ncov){ # THIS IS for aftplot!!
    ##    uppe <- exp(-sum(new.data[1:ncov] * x$coefficients[1:ncov]) / p)
    ##    lambda <- lambda * uppe
    ##}
    if (is.null(xlim)){
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))
        if (xlim[1] <= 0) xlim[1] <- 0.001 * xlim[2] ## Avoid Inf at 0! (hack!?)
    }
    npts <- 4999
    xx <- seq(xlim[1], xlim[2], length = npts)
    haz <- matrix(0, ncol = npts, nrow = ns)
    sur <- haz
    Haz <- haz
    ## hazard
    if (x$dist == "weibull"){
        for (i in 1:ns){
            scal <- scale[i] / (score)^(1/shape[i])
            haz[i, ] <- hweibull(xx, shape = shape[i],
                                 scale = scal)
            sur[i, ] <- pweibull(xx, shape = shape[i],
                                 scale = scal, lower.tail = FALSE)
            Haz[i, ] <- Hweibull(xx, shape = shape[i],
                                 scale = scal)
        }
        dist <- "Weibull"
    }else if (x$dist == "pch"){
        for (i in 1:ns){
            haz[i, ] <- hpch(xx, x$cuts, score * x$hazards[i, ])
            sur[i, ] <- ppch(xx, x$cuts, score * x$hazards[i, ],
                             lower.tail = FALSE)
            Haz[i, ] <- Hpch(xx, x$cuts, score * x$hazards[i, ])
        }
        dist <- "Pcwise const"
    }else if (x$dist == "loglogistic"){
        for (i in 1:ns){
            haz[i, ] <- hllogis(xx, shape = shape[i],
                                 scale = scale[i]) * score
            sur[i, ] <- pllogis(xx, shape = shape[i],
                                 scale = scale[i], lower.tail = FALSE)^score
            Haz[i, ] <- Hllogis(xx, shape = shape[i],
                                 scale = scale[i]) * score
        }
        dist <- "Loglogistic"

    }else if (x$dist == "lognormal"){
        sdlog <- 1 / shape
        meanlog <- log(scale)
        for (i in 1:ns){
            haz[i, ] <- hlnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i]) * score
            sur[i, ] <- plnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i], lower.tail = FALSE)^score
            Haz[i, ] <- Hlnorm(xx, meanlog = meanlog[i],
                               sdlog = sdlog[i]) * score
        }
        dist = "Lognormal"
    }else if (x$dist == "ev"){
        for (i in 1:ns){
            haz[i, ] <- hEV(xx, shape = shape[i],
                                 scale = scale[i]) * score
            sur[i, ] <- pEV(xx, shape = shape[i],
                                 scale = scale[i], lower.tail = FALSE)^score
            Haz[i, ] <- HEV(xx, shape = shape[i],
                                 scale = scale[i]) * score
        }

        dist = "Extreme value"
    }else if (x$dist == "gompertz"){
        if (x$param == "canonical"){
            for (i in 1:ns){
                haz[i, ] <- hgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      param = "canonical")## * score
                sur[i, ] <- pgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      lower.tail = FALSE,
                                      param = "canonical") ##^score
                Haz[i, ] <- Hgompertz(xx, shape = score * shape[i],
                                      scale = scale[i],
                                      param = "canonical") ##* score
            }
        }else if (x$param == "rate"){
            ##scale <- log(scale) # Added 15 Jul 2018. Removed 19 Jul 2018!
            ##scale <- 1 / scale
            scale <- log(scale)
            ##if (printlevel){
              ##  cat("scale = ", scale, "\n")
                ##cat("shape = ", shape, "\n")
            ##}
            for (i in 1:ns){
                ##haz[i, ] <- hgompertz(xx, shape = score * shape[i], scale = scale[i])
                haz[i, ] <- shape[i] * exp(xx * scale[i]) * score
                ##Haz[i, ] <- Hgompertz(xx, shape = score * shape[i], scale = scale[i])
                Haz[i, ] <- shape[i] * score * expm1(xx * scale[i]) / scale[i]
                sur[i, ] <- exp(-Haz[i, ]) # Fixed 15 Jul 2018 (order [i, ])
            }
        }

        dist = "Gompertz"
    }

    if ("haz" %in% fn){

        if (is.null(ylim)) {
            ylim0 <- c(0, max(haz))
        }else{
            ylim0 <- ylim
        }
        ##if (min(p) < 1) ylim0[2] <- min(ylim0[2], max(haz[, -1]))

        ##if (is.null(xlab)) xlab <- "Duration"
        if (is.null(ylab)) ylab <- "Hazard"
        if (is.null(main)){
            hmain <- paste(dist, "hazard function")
        }else{
            hmain <- main
        }
        plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             col = col[1], lty = lty[1],
             xlab = xlab, ylab = ylab, main = hmain, ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, haz[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
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
    }
    ## Cumulative hazard
    if ("cum" %in% fn && fig){
        
        if (is.null(ylim)){
            ylim0 <- c(0, max(Haz))
        }else{
            ylim0 <- ylim
            ylim0[2] <- max(ylim0[2], max(Haz))
        }
        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Cumulative Hazards"
        if (is.null(main)){
            Hmain <- paste(dist, "cumulative hazard function")
        }else{
            Hmain <- main
        }
        plot(xx, Haz[1, ], type = "l", xlim = xlim, ylim = ylim0,
             xlab = xlab, ylab = ylab, main = Hmain, col = col[1],
             lty = lty[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, Haz[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        ##abline(h = 0)
        ##abline(v = 0)
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = "topleft",  legend = x$strata, lty = lty,
                       col = col, inset = 0.001)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend, legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
        
    }

    ## density
    if ("den" %in% fn){

        den <- haz * sur
        ##if (is.null(ylim))
        ylim <- c(0, max(den))

        ##if (min(p) < 1) ylim[2] <- min(max(den[, -1]))

        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Density"
        if (is.null(main)){
            dmain <- paste(dist, "density function")
        }else{
            dmain <- main
        }
        plot(xx, den[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = dmain, lty = lty[1],
             col = col[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, den[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        abline(h = 0)
        abline(v = 0)
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = printLegend,  legend = x$strata, lty = lty,
                       inset = 0.001,
                       col = col)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend,  legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
        
    }
    
    ## Survivor function
    if ("sur" %in% fn){
        
        ##if (is.null(ylim))
        ylim <- c(0, 1)

        ##if (is.null(xlab))
        ##xlab <- "Duration"
        ##if (is.null(ylab))
        ylab <- "Survival"
        if (is.null(main)){
            smain <- paste(dist, "survivor function")
        }else{
            smain <- main
        }
        plot(xx, sur[1, ], type = "l", xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = smain, lty = lty[1],
             col = col[1], ...)
        if (ns > 1){
            for (i in 2:ns){
                lines(xx, sur[i, ], type = "l", lty = lty[i],
                      col = col[i], ...) # ', ...' added in 2.4-4
            }
        }
        abline(h = 0)
        abline(v = 0)
        if (is.logical(printLegend)){
            if ((ns > 1) && printLegend){
                legend(x = "bottomleft",  legend = x$strata, lty = lty,
                       inset = 0.001,
                       col = col)
            }
        }else{
            if ((ns > 1) && is.character(printLegend)){
                    legend(x = printLegend,  legend = x$strata,
                           lty = lty, inset = 0.001, col = col)
            }
        }
        
    }

    if (fn == "cum" & !fig) invisible(list(xx = xx, Haz = Haz))
    ##par(oldpar)
}
