plot.weibreg <- function(x,
                         main = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         new.data = x$means,
                         ...){
    if (!inherits(x, "weibreg")) stop("Works only with 'weibreg' objects.")
    if (x$pfixed) stop("True exponential hazards are not plotted")
    ncov <- length(x$means)
    ns <- x$n.strata
    lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    p <- exp(x$coefficients[ncov + (1:ns) * 2])
    if (is.null(xlim))
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))
    
    npts <- 199
    xx <- seq(xlim[1], xlim[2], length = npts)
    if (xx[1] <= 0) xx[1] <- 0.001
    haz <- matrix(ncol = npts, nrow = ns)

    for (i in 1:ns){
        tl <- xx / lambda[i]
        haz[i, ] <- (p[i] / lambda[i]) *
            tl^(p[i]-1) * exp(new.data[1:ncov] * x$coefficients[1:ncov])
    }

    if (is.null(ylim)) ylim <- c(0, max(haz))
    if (is.null(xlab)) xlab <- "Duration"
    if (is.null(ylab)) ylab <- "Hazard"
    if (is.null(main)) main <- "Weibull hazard function"
    plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ...)
    if (ns > 1){
        for (i in 2:ns){
            lines(xx, haz[i, ], type = "l", lty = i)
        }
    }
    abline(h = 0)
}
