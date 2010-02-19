check.dist <- function(sp, pp, new.data = sp$means, main = NULL){
    if (!inherits(sp, "coxreg"))
        stop ("First argument must be of type 'coxreg'")
    if (!inherits(pp, "phreg")) stop ("Second argument must be of type 'phreg'")

    x <- plot.coxreg(sp, fn = "cum", fig = FALSE, new.data = new.data)$x
    y.max <- max(x[[1]][, 2])
    if (length(x) > 1){
        for (i in 2:length(xx)) y.max <- max(y.max, x[[i]][, 2])
    }        
    plot.phreg(pp, fn = "cum", fig = TRUE, new.data = new.data,
               ylim = c(0, y.max), main = main)
    for (rr in 1:length(x)){
        xx <- x[[rr]]
        lines(xx[, 1], xx[, 2], type = "s", lty = 2, col = "red")
    }
}

