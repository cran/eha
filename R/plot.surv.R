plot.Surv <- function(x,  ## A survival object
                      strata = rep( 1, length(exit) ),
                      limits = FALSE,
                      conf = 0.95,
                      xlim = NULL,
                      ylim = NULL,
                      main = "Survivor function(s)",
                      xlab = "Duration",
                      ylab = "Remaining fraction",
                      ...)
  {
    require(survival)
    ## Input data:
    ##
    ## enter : left truncation point
    ## exit  : exit time point
    ## event : if zero, a censored observation; otherwise an event.
    ## strata : one curve for each value of strata.
    ## limits: if TRUE, and only one strata, pointwise confidence
    ##         limits (Greenwoods formula, log(-log) type.
    ## conf  : confidence level. Can be given as a percentage.
    
    ## Check input data:
    if (!inherits(x, "Surv")) 
      stop("First arg must be of type 'Surv'")

    if (ncol(x) == 3){
      enter <- x[, 1]
      exit <- x[, 2]
      event <- x[, 3]
      n <- length(exit)
    }else{
      exit <- x[, 1]
      n <- length(exit)
      enter <- rep(0, n)
      event <- x[, 2]
    }

    if (is.na(strata)) strata <- rep(1, n)
    if (length(enter) != n)stop("enter and exit must have equal length.")
    if (length(event) != n)
      stop("event and exit must have equal length.")
    if (length(strata) != n)
      stop("strata and exit must have equal length.")
    if (min(exit - enter) <= 0) stop("Interval length must be positive.")
    if (conf > 1) conf <- conf / 100 ## conf given as a percentage(?)
    if ( (conf < 0.5) | (conf >=1) ) stop("Bad conf value")

    grupp <- as.character(strata)
   
    strata <- sort(unique(grupp))
    no.of.groups <- length(strata)
    if (no.of.groups > 9)
      stop("Too many groups. No more than 9 are allowed.")
    
    ##
    if (length(strata) > 1) limits <- FALSE # No limits if multiple curves.

    ## Check xmin, xmax:

    if ( is.null(xlim) ){
      x.max <- max(exit) ## Must be better?!
      x.min <- 0
      xlim = c(0, x.max)
    }else{
      x.min <- xlim[1]
      x.max <- xlim[2]
    }

    
    if ( is.null(ylim) ){
      ylim = c(0, 1)
    }
      
    gang <- 0

    for (stratum in strata)
      {
        atom <- table.events(enter[grupp == stratum],
                             exit[grupp == stratum],
                             event[grupp == stratum])
        
        gang <- gang + 1
        surv <- c( 1, cumprod(1 - atom$events / atom$riskset.sizes) )
        if (gang == 1)
          {
            X <- rep(c(0, atom$times), each = 2)[-1]
            Y <- rep(surv, each = 2)[-2*length(surv)]
            plot(X, Y, type = "l",
                 xlab = xlab, ylab = ylab,
                 main = main, xlim = xlim, ylim = ylim,
                 lty = gang%%no.of.groups + 1, ...)
            if (limits) ## Greenwood's formula,
                        ## Kalbfleisch & Prentice, p. 15 (note error!).
              {
                q.alpha <- abs(qnorm((1 - conf) / 2))
                survived <- (atom$riskset.size - atom$events)
                se <- sqrt(cumsum(atom$events /
                                  ( atom$riskset.sizes * survived )
                                  )
                           )/
                            cumsum(-log(survived / atom$riskset.sizes))
                upper <- surv ^ exp(q.alpha * c(0, se))
                lower <- surv ^ exp(-q.alpha * c(0, se))
                X <- rep(c(0, atom$times), each = 2)[-1]
                Y <- rep(upper, each = 2)[-2*length(upper)]

                lines(X, Y, type = "l",
                      lty = gang%%no.of.groups + 2)
                Y <- rep(lower, each = 2)[-2*length(lower)]
                lines(X, Y, type = "l",
                      lty = gang%%no.of.groups + 2)
              }
          }
        else
          {
            X <- rep(c(0, atom$times), each = 2)[-1]
            Y <- rep(surv, each = 2)[-2*length(surv)]
            lines(X, Y, type = "l", 
                  lty = gang%%no.of.groups + 1)
          }
      }
    abline(h = 0)
    abline(v = 0)
    if (no.of.groups > 1)
      {
        colors <- (1:no.of.groups)%%no.of.groups + 1
        legend(x.min, 0, xjust = 0, yjust = 0,
               legend = strata, lty = colors)
      }
  }
