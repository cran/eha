check.dist <- function(sp, pp,
                       new.data = sp$means){
    if (!inherits(sp, "coxreg"))
        stop ("First argument must be of type 'coxreg'")
    if (!inherits(pp, "phreg")) stop ("Second argument must be of type 'phreg'")

    plot.phreg(pp, fn = "cum", fig = TRUE, new.data = new.data)
    x <- plot.coxreg(sp, fn = "cum", fig = FALSE, new.data = new.data)$x[[1]]
    lines(x[, 1], x[, 2], type = "s", lty = 2, col = "red")
}

