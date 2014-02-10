plot.coxreg <- function(x,
                        fn = c("cum", "surv", "log", "loglog"),
                        fig = TRUE,
                        xlim = NULL,
                        ylim = NULL,
                        main = NULL,
                        xlab = "Duration",
                        ylab = "",
                        col,  ## New 5 dec 2013.
                        lty,  ## New 17 Jan 2014
                        printLegend = TRUE, ## New 17 Jan 2014
                        new.data = NULL,
                        ...){

    if (missing(col)) col <- "black" ## New 2013-12-05
    if (missing(lty)) lty <- 1:length(x$hazards) # No. of strata
    if (!is.null(new.data)) warning("argument 'newdata' is not used any more")
    fn <- fn[1]
    if (!inherits(x, "coxreg")) stop("Works only with 'coxreg' objects.")
    if (is.null(x$hazards)){
        cat("No 'hazards' object found. Must be fixed!!!!")
        return(NULL)
    }
    if (!(fn %in% c("cum", "surv", "log", "loglog")))
        stop(paste(fn, "is an illegal value of 'fn'"))

    if (FALSE){ # From 2.4-0: hazards are hazards! (But I may change my mind)
    ##if ((!is.null(new.data)) && (!is.null(x$coefficients))){
        score <- exp(sum((new.data - x$means) * x$coefficients))
        if (is.data.frame(x$hazards)){
            x$hazards$hazard <- x$hazards$hazard * score
        }else{
            for (i in 1:length(x$hazards))
                x$hazards[[i]][, 2] <- 1 - (1 - x$hazards[[i]][, 2])^score
            ## Depends on 'hazards.f' K&p, II, p. 116. This must be
            ## reconsidered in the future!!
            ##x$hazards[[i]][, 2] <- score * x$hazards[[i]][, 2]
        }
    } # end{ if (FALSE) }
    if (is.data.frame(x$hazards)){
        plot.hazards(x$hazards, fn, fig, xlim, ylim, main, xlab, ylab, ...)
    }else{
        plot.hazdata(x$hazards, x$strata,
                     fn, fig, xlim, ylim, main, xlab, ylab,
                     col = col, lty = lty, printLegend = printLegend, ...)
    }
}

