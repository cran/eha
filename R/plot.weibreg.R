## Doesn't work with 'newdata'! Will be fixed.

plot.weibreg <- function(x, new.data = rep(0, length(x$means)), ...){
    if (!inherits(x, "weibreg")) stop("Works only with 'weibreg' objects.")
    ncov <- length(x$means)
    ns <- x$n.strata
    lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    p <- exp(x$coefficients[ncov + (1:ns) * 2])
    xlim <- c(min(x$y[, 1]), max(x$y[, 2]))
    
    npts <- 199
    xx <- seq(xlim[1], xlim[2], length = npts)
    if (xx[1] <= 0) xx[1] <- 0.001
    haz <- matrix(ncol = npts, nrow = ns)

    for (i in 1:ns){
        tl <- xx / lambda[i]
        haz[i, ] <- (p[i] / lambda[i]) *
            tl^(p[i]-1)## * exp(new.data[1:ncov] * x$coefficients[1:ncov])
    }

    ylim <- c(0, max(haz))

    plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim,
         xlab = "age", ylab = "hazard", main = "Webull hazards")
    if (ns > 1){
        for (i in 2:ns){
            lines(xx, haz[i, ], type = "l", lty = i)
        }
    }
}
