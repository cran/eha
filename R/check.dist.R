check.dist <- function(sp, pp, main = NULL){
    if (!inherits(sp, "coxreg"))
        stop ("First argument must be of type 'coxreg'")
    if (!inherits(pp, "phreg"))
        stop ("Second argument must be of type 'phreg' or 'pchreg'")

    x.max <- max(pp$y[, 2])
    x <- plot.coxreg(sp, fn = "cum", fig = FALSE)$x
    if (is.null(x)){
        cat("must be fixed in check.dist! or plot.coxreg")
        return(x)
    }
    y.max <- max(x[[1]][, 2])
    if (length(x) > 1){
        for (i in 2:length(x)) y.max <- max(y.max, x[[i]][, 2])
    }
    plot(pp, fn = "cum", fig = TRUE,
         ylim = c(0, y.max), main = main)
    for (rr in 1:length(x)){
        xx <- x[[rr]]
        xx <- rbind(xx, c(x.max, xx[NROW(xx), 2])) # Added 2011-08-10 (2.0-3)
        lines(xx[, 1], xx[, 2], type = "s", lty = 2, col = "red")
    }
}

