plot.coxreg <- function(x,
                        fn = c("cum", "surv", "log", "loglog"),
                        fig = TRUE,
                        xlim = NULL,
                        ylim = NULL,
                        main = NULL,
                        xlab = "Duration",
                        ylab = "",
                        new.data = x$means,
                        ...){
    fn <- fn[1]
    if (!inherits(x, "coxreg")) stop("Works only with 'coxreg' objects.")
    if (is.null(x$hazards)) stop("No 'hazards' object found.")
    if (!(fn %in% c("cum", "surv", "log", "loglog")))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if ((!is.null(new.data)) && (!is.null(x$coefficients))){
        score <- exp(sum((new.data - x$means) * x$coefficients))
        for (i in 1:length(x$hazards))
            x$hazards[[i]][, 2] <- 1 - (1 - x$hazards[[i]][, 2])^score
            ## Depends on 'hazards.f' K&p, II, p. 116. This must be
            ## reconsidered in the future!!
            ##x$hazards[[i]][, 2] <- score * x$hazards[[i]][, 2]
    }
    plot.hazdata(x$hazards, fn, fig, main, xlim, ylim, xlab, ylab, ...)

}

