plot.cum <- function(x, ## A survival object
                     group = rep(1, length(exit)),
                     main = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     log.scale = FALSE,
                     ...
                     )
  {
      warning("'plot.cum' is deprecated (use 'plot.Surv')")
      if (log.scale) foo <- "loglog"
      else foo <- "cum"

      plot.Surv(x, main = main, xlab = xlab,
                ylab = ylab, fn = foo, ...)
  }
