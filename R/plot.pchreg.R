plot.pchreg <- function(x,
                        fn = c("haz", "cum", "den", "sur"),
                        main = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        xlab = "Duration",
                        ylab = NULL,
                        new.data = NULL,
                         ...){

    if (is.null(new.data)) new.data <- x$means
    
    if (!inherits(x, "pchreg")) stop("Works only with 'pchreg' objects.")
    ##if (x$pfixed) stop("True exponential hazards are not plotted")
    if (!(all(fn %in% c("haz", "cum", "den", "sur"))))
        stop(paste(fn, "is an illegal value of 'fn'"))


    if (length(fn) >= 3){
        oldpar <- par(mfrow = c(2, 2))
        on.exit(par(oldpar))
    }else if (length(fn) == 2){
        oldpar <- par(mfrow = c(2, 1))
        on.exit(par(oldpar))
    }
    ncov <- length(x$means)

    if (ncov){
        score <- exp(sum((new.data - x$means) * x$coefficients[1:ncov]))
    }else{
        score <- 1
    }

    ##if (ncov){ # THIS IS for aftplot!!
    ##    uppe <- exp(-sum(new.data[1:ncov] * x$coefficients[1:ncov]) / p)
    ##    lambda <- lambda * uppe
    ##}
    if (is.null(xlim))
        xlim <- c(min(x$y[, 1]), max(x$y[, 2]))

    npts <- 4999
    xx <- seq(xlim[1], xlim[2], length = npts)
    ##if (xx[1] <= 0) xx[1] <- 0.001
    yy <- numeric(length(xx))
    cuts <- c(xlim[1], x$cuts, xlim[2])
    haz <- x$hazards
    n.ivl <- length(cuts) - 1
    for (i in 1:n.ivl){
        yy[(xx > cuts[i]) & (xx <= cuts[i+1])] <- x$hazards[i]
    }
    skal <- NULL
    ## hazard

    dist <- "Pch"
    if ("haz" %in% fn){

        if (is.null(ylim)){
            hylim <- c(0, max(yy))
        }else{
            hylim <- ylim
        }

        if (is.null(xlab)){
            hxlab <- "Duration"
        }else{
            hxlab <- xlab
        }
        if (is.null(ylab)){
            hylab <- "Hazards"
        }else{
            hylab = ylab
        }
        if (is.null(main)){
            hmain <- paste(dist, "hazard function")
        }else{
            hmain <- main
        }
        plot(cuts[1:2], c(haz[1], haz[1]), type = "l",
             xlim = xlim, ylim = hylim,
             xlab = hxlab, ylab = hylab, main = hmain, ...)
        for (i in 2:(length(cuts) - 1)){
            lines(cuts[i:(i+1)], c(haz[i], haz[i]))
        }
             
        ##plot(xx, yy, ##type = "l",
          ##   xlim = xlim, ylim = ylim,
        abline(h = 0)
        abline(v = 0)
    }
    ## Cumulative hazard
    if ("cum" %in% fn){

        Haz <- numeric(length(cuts))
        Haz[1] <- 0
        for (i in 2:length(Haz)){
            Haz[i] <- Haz[i-1] + haz[i-1] * (cuts[i] - cuts[i-1])
        }
        if (is.null(ylim)){
            hylim <- c(0, max(Haz))
        }else{
            hylim <- ylim
            hylim[2] <- max(hylim[2], max(Haz))
        }
        if (is.null(xlab)){
            hxlab <- "Duration"
        }else{
            hxlab <- xlab
        }
        if (is.null(ylab)){
            hylab <- "Cumulative Hazards"
        }else{
            hylab <- ylab
        }
        if (is.null(main)){
            hmain <- paste(dist, "cumulative hazards function")
        }else{
            hmain <- main
        }
        plot(cuts, Haz, type = "l",
             xlim = xlim, ylim = hylim,
             xlab = hxlab, ylab = hylab, main = hmain, ...)

        abline(h = 0)
        abline(v = 0)
    }
    ## density
    if ("den" %in% fn){
        
        ##if (is.null(ylim))
        yy <- dpch(xx, x$cuts, x$hazards)
        if (is.null(xlab)){
            hxlab <- "Duration"
        }else{
            hxlab <- xlab
        }
        if (is.null(ylab)){
            hylab <- "Density"
        }else{
            hylab <- ylab
        }
        if (is.null(main)){
            hmain <- paste(dist, "density function")
        }else{
            hmain <- main
        }
        if (is.null(ylim)){
            hylim <- c(0, max(yy))
        }else{
            hylim <- ylim
        }
        who <- (xx >= 0) & (xx <= x$cuts[1]) 
        plot(xx[who], yy[who], type = "l", xlim = xlim, ylim = hylim,
             xlab = hxlab, ylab = hylab, main = hmain, ...)
        for (i in 2:(length(cuts) - 1)){
            who <- (xx > cuts[i]) & (xx <= cuts[i + 1]) 
            lines(xx[who], yy[who])
        }
        
        abline(h = 0)
        abline(v = 0)
    }
    ## Survivor function
    if ("sur" %in% fn){

        yy <- ppch(xx, x$cuts, x$hazards, lower.tail = FALSE)

        ##if (is.null(ylim))
        ylim <- c(0, 1)

        if (is.null(xlab)){
            hxlab <- "Duration"
        }else{
            hxlab <- xlab
        }
        if (is.null(ylab)){
            hylab <- "Survival"
        }else{
            hylab <- ylab
        }
        if (is.null(main)){
            hmain <- paste(dist, "survivor function")
        }else{
            hmain <- main
        }
        plot(xx, yy, type = "l", xlim = xlim, ylim = ylim,
             xlab = hxlab, ylab = hylab, main = hmain, ...)
        abline(h = 0)
        abline(v = 0)
    }
    ##par(oldpar)
}
