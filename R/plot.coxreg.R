#' Plot method for \code{coxreg} objects
#' 
#' A plot of a baseline function of a \code{coxreg} fit is produced, one curve
#' for each stratum. A wrapper for \code{\link[survival]{plot.survfit}}.
#' 
#' @param x A \code{coxreg} object
#' @param fn What should be plotted? Default is "cumhaz", and the other choices
#' are "surv", "log", and "loglog".
#' @param conf.int logical or a value like 0.95 (default for one curve).
#' @param fig logical. If \code{TRUE} the plot is actually drawn, otherwise
#' only the coordinates of the curve(s) are returned.
#' @param xlim Start and end of the x axis.
#' @param ylim Start and end of the y axis.
#' @param main A headline for the plot
#' @param xlab Label on the x axis.
#' @param ylab Label on the y axis.
#' @param col Color of the curves. Defaults to 'black'.
#' @param lty Line type(s).
#' @param printLegend Either a logical or a text string; if \code{TRUE}, a
#' legend is printed at a default place, if \code{FALSE}, no legend is printed.
#' Otherwise, if a text string, it should be one of "bottomleft",
#' "bottomright", "topleft", etc., see \code{\link{legend}} for all possible
#' choices.
#' @param ... Other parameters to pass to the plot.
#' @return An object of class \code{hazdata} containing the coordinates of the
#' curve(s).
#' @export
plot.coxreg <- function(x,
                        fn = c("cum", "surv", "log", "loglog"),
                        conf.int = FALSE,
                        fig = TRUE,
                        xlim = NULL,
                        ylim = NULL,
                        main = NULL,
                        xlab = "Duration",
                        ylab = "",
                        col = 1,
                        lty = 1, 
                        printLegend = TRUE,
                        ...){
   if (!inherits(x, c("coxreg"))){
      stop("Works only with 'coxreg' objects")
   }
   ##if (x$method %in% c("ml", "mppl")){
   if (is.null(x$hazards)){
      x$hazards <- getHaz(x$y, strats = x$stratum, 
                          score = exp(x$linear.predictors))
   }
   ##if (!is.null(x$hazards)){
   if (all(names(x$hazards) == as.character(1:length(x$hazards)))){
      if (!is.null(x$strata)){
         names(x$hazards) <- x$strata 
      }
   }
   op <- par(las = 1)
   if (!conf.int){
      res <- plot.hazdata(x$hazards, fn = fn, fig = fig, xlim = xlim, ylim = ylim, 
                     main = main, xlab = xlab, ylab = ylab, col = col, 
                      lty = lty, printLegend = printLegend)
      par(op)
      return(invisible(res))
   }
   if (conf.int){ # Old stuff kept as backup...
      x$means <- numeric(length(x$coefficients))
      class(x) <- "coxph"
      fn <- fn[1]
      if (fn == "cum"){
         fn = "cumhaz"
      }else{
         if (fn == "loglog"){
            fn = "cloglog"
         }
      }
      
      if (is.null(xlim)) {
         if (NCOL(x$y) == 3){
            xlim <- c(min(as.numeric(x$y[, 1])), max(as.numeric(x$y[, 2])))
         }
      }
      if (is.null(ylim)){
         
      }
      
      if (FALSE){
         ##if (!is.null(x$nullModel) && x$nullModel){ # NOT NEEDED!!!
         ## Null model and from coxreg...
         plot(x$y ~ 1, fun = fn, xlab = xlab, ylab = ylab, main = main, 
              xlim = xlim, col = col, lty = lty, conf.int = conf.int)
         if (printLegend && length(x$stratum)){
            legend("bottom", legend = x$stratum, lty = lty, col = col)
         }
         
      }else{
         res <- survival::survfit(x)
         plot(res, fun = fn, xlab = xlab,
              ylab = ylab, main = main, xlim = xlim, col = col, lty = lty,
              conf.int = conf.int)
         if (printLegend && length(x$stratum) >= 2){
            legend("bottom", legend = x$stratum, lty = lty, col = col)
         } 
      }
   } # End if (FALSE) "old stuff"
   invisible(res)
}
