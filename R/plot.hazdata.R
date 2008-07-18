plot.hazdata <- function(x,
                         fn = c("cum", "surv", "log", "loglog"),
                         fig = TRUE,
                         xlim = NULL,
                         ylim = NULL,
                         main = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         ...
                         ){
    ## x is of type 'hazdata', which is a list with two-column matrices
    ## as components, one component per stratum. The first column contains
    ## risktimes, and the second column the corresponding 'hazard atoms'.

    if (!inherits(x, "hazdata")){
        if (!inherits(x, "coxreg")){
            stop("First argument must be of type 'hazdata' or 'coxreg'")
        }else{
            y <- x
            x <- x$hazards
            if (is.null(x)) stop("No 'hazards' component present")
        }
    }
    if (!(fn %in% c("cum", "surv", "log", "loglog")))
        stop(paste(fn, "is an illegal value of 'fn'"))

    n.strata <- length(x)
    fn <- fn[1]

    yVal <- function(x){
        if (fn == "cum") return(cumsum(x))
        if (fn %in% c("log", "loglog")) return(log(cumsum(x)))
        n <- length(x)
        s <- numeric(n)
        s[1] <- 1 - x[1]
        if (n > 1){
            for (rs in 2:n){
                s[rs] <- s[rs - 1] *(1 - x[rs])
            }
        }
        return(s)
    }

    max.x <- max(x[[1]][, 1])
    min.x <- min(x[[1]][, 1])

    if (fn == "loglog"){
        max.x <- log(max.x)
        min.x <- log(min.x)
    }

    max.y <- -1000 # What else :-)
    min.y <- 0
    for (i in 1:n.strata){
        x[[i]][, 2] <- yVal(x[[i]][, 2])
        max.y <- max(c(max.y, x[[i]][, 2]))
        min.y <- min(c(min.y, x[[i]][, 2]))
        if (fn == "loglog") x[[i]][, 1] <- log(x[[i]][, 1])
        max.x <- max(c(max.x, x[[i]][, 1]))
        min.x <- min(c(min.x, x[[i]][, 1]))
    }

    if (is.null(xlim)) {
        xlim <- c(min.x, max.x)
        if (fn != "loglog") xlim[1] <- 0.0
    }
    if (length(xlim) != 2) stop("'xlim' must be a vector of length two")

    if (is.null(ylim)){
        ylim <- c(min.y, max.y)
        if ( (fn == "surv") || (fn == "cum") ) ylim[1] <- 0.0
    }
    if (length(ylim) != 2) stop("'ylim' must be a vector of length two")

    if (fn != "loglog"){
        for (i in 1:n.strata){
            if (fn == "surv"){
                x[[i]] <- rbind(c(x[[i]][1, 1], 1), x[[i]])
            }else{
                x[[i]] <- rbind(c(x[[i]][1, 1], 0), x[[i]])
            }
            ##x[[i]][, 1] <- c(x[[i]][1, 1], x[[i]][, 1])
        }
    }
    if (fig){
        plot(x[[1]][, 1], x[[1]][, 2], type = "s",
             xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, main = main, lty = 1, ...)
        if (n.strata > 1){
            for (i in 2:n.strata){
                lines(x[[i]][, 1], x[[i]][, 2], type = "s", lty = i)
            }
        }
        abline(h = 0)
    }
    invisible(list(x = x, fn = fn))
}

