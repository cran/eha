plot.cum <- function(x, ## A survival object
                     group = rep(1, length(exit)),
                     main = "Cumulative hazards functions",
                     xlab = "Duration",
                     ylab = "",
                     log.scale = FALSE, ...
                     )
  {

    require(survival)
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

 
    if (length(enter) != n) stop("enter and exit must have equal length.")
    if (length(event) != n) stop("event and exit must have equal length.")
    if (length(group) != n) stop("group and exit must have equal length.")
    if (min(exit - enter) <= 0) stop("Interval lengths must be positive.") 

    if (is.factor(group)){
      strata <- levels(group)
    }else{
      group <- as.character(group)
      strata <- sort(unique(group))
    }
    no.of.groups <- length(strata)
    if (no.of.groups > 9)
      stop("Too many groups. No more than 9 are allowed.")
     
    ## Check for ylim, xlim in coming plots:
    y.max <- 0
    y.min <- 1
    x.max <- 0
    x.min <- 1e103
    for (stratum in strata)
      {
        atom <- table.events(enter[group == stratum],
                        exit[group == stratum],
                        event[group == stratum])
        y.max <- max( y.max, sum(atom$events / atom$riskset.sizes) )
        y.min <- min(y.min, atom$events[1] / atom$riskset.sizes[1])
        x.max <- max(x.max, atom$times)
        x.min <- min(x.min, atom$times[1])
      }

    ## Start plotting:        
    gang <- 0

    for (stratum in strata)
      {
        atom <- table.events(enter[group == stratum],
                             exit[group == stratum],
                             event[group == stratum])
        
        gang <- gang + 1

        if (log.scale)
          {
            xy <- "xy"
            if (x.min <= 0) x.min <- x.max / 100
            cum <- c(cumsum(atom$events / atom$riskset.sizes))
          }
        else
          {
            xy <- ""
            y.min <- 0
            atom$times <- c(0, atom$times) 
            cum <- c(0, cumsum(atom$events / atom$riskset.sizes))
          }
        n.po <- length(cum)
        x.po <- c(atom$times[1], rep(atom$times[-1], rep(2, n.po - 1)))
        y.po <- c(rep(cum[-n.po], rep(2, n.po - 1)), cum[n.po])
        if (gang == 1)
          {
            plot(x.po, y.po, type = "l",
                 xlab = xlab, ylab = ylab, 
                 main = main, log = xy,
                 xlim = c(x.min, x.max), ylim = c(y.min, y.max),
                 lty = gang%%(no.of.groups + 1), ...)
            abline(h = 0)
          }
        else
          {
            lines(x.po, y.po, type = "l",  
                  lty = gang%%(no.of.groups + 1))
          }
        if (no.of.groups > 1)
          {
            colors <- (1:no.of.groups)%%(no.of.groups + 1)
            legend(x.min, y.max, lty = colors, legend = strata)
          }
      }
  }
