# Puts in 'communal time dependent covariates by 'splitting spells'.
# Göran Broström (2001).

make.communal <- function(dat,
                          com.dat,
                          ##com.info,
                          communal = TRUE,
                          start,
                          period = 1,
                          lag = 0,
                          surv = c("enter", "exit", "event", "birthdate"),
                          tol = 0.0001,
                          fortran = TRUE){
  
  ## 'dat' is a data frame with variables:
  ## birthdate = birth date
  ## enter = left truncation time
  ## exit  = right censoring/event time
  ## event = event indicator (0 if no event).
  ## other covariates.

  ## 'com.dat' is a data frame with columns communal covariates.

  ## Formula: datum = bdate + risktime * as.numeric(communal) + lag
  ## Note: 'scale' may only be 0 or 1. If scale = 1,
  ## the 'lag' must be <= 0! ("causality")
  ## Leaving:
  ## either 'datum = bdate + risktime + lag' (communal == TRUE, lag == 0)
  ## or     'datum = bdate + lag'      (communal == FALSE).

  ## NOTE: names(com.dat) must be != names(dat) !!!

  if (!is.data.frame(dat))stop("dat must be a data frame")
  if (!is.data.frame(com.dat))stop("com.dat must be a data frame")
  ##if (!is.vector(com.info))stop("com.info must be a vector")
  ##if (length(com.info) != 4)stop("com.info must be a vector of length 4")
  ##scale <- com.info[4]
  ##if ( !((scale == 0) || (scale == 1)) ) stop("scale must be 0 or 1")
  ## lag <- com.info[3]
  
  if (length(surv) != 4) stop("surv must have length 4")
  fixed.names <- names(dat)
  surv.indices <- match(surv, fixed.names)
  if (length(which(is.na(surv.indices)))){
    x <- which(is.na(surv.indices))
    stop(paste(surv[x], " is not a name in the fixed data frame."))
  }
  com.names <- names(com.dat)
  if ( length(x <- which(!is.na(match(com.names, c(surv, fixed.names))))) )
    stop(paste(com.names[x], "are names in fixed data frame.")) 
  
  nn <- nrow(dat)
  n.years <- nrow(com.dat)
  n.com <- NCOL(com.dat) ## No. of communal covariates.  
  ## Function to calculate start and stop period

  if (n.com > 1) stop("Only one communal covariate at a time!")
  ## NOTE: Will only work with n.com == 1 for now!!

  cuts <- start + c(0, (1:n.years) * period) - lag

  beg.per <- cuts[1]
  end.per <- cuts[n.years + 1]
  iv.length <- period
  

  ## cut off in calendar time:
  spell.tot <- sum(dat[, surv.indices[2]] - dat[, surv.indices[1]])
  dat <- cal.window(dat, c(beg.per, end.per), surv)
  if (sum(dat[, surv.indices[2]] - dat[, surv.indices[1]]) < spell.tot)
    warning("Spells are cut")
  nn <- nrow(dat)
    
  if (!communal){ ## "Fixed" communal!

    get.per <- function(dates)
      pmin(pmax(1, ceiling((dates - beg.per) / iv.length)),
           n.years)
    dates <- dat[, surv.indices[4]]
    ppp <- get.per(dates)
    yy <- matrix(0, ncol = n.com, nrow = nn)
    for (i in 1:n.com){
      yy[, i] <- com.dat[ppp, i]
    }
    yy <- as.data.frame(yy)
    names(yy) <- com.names
    yy <- cbind(dat, yy)
    
  }else{ ## Real communal!

    get.iv <- function(dates)

      cbind(pmin(pmax( 1, floor((dates[, 1, drop = FALSE] - beg.per) /
                          iv.length) + 1 ), n.years ),
            pmin(pmax( 1, ceiling((dates[, 2, drop = FALSE] - beg.per) /
                                 iv.length) ), n.years) 
	)

    event.boolean <- is.logical(dat[, surv.indices[3]])
    xx <- cbind(dat[, surv.indices, drop = FALSE], 1:nn)
    xx[, 3] <- as.numeric(xx[, 3]) ## Could be a boolean vector.
    xx <- as.matrix(xx)
    if (!is.numeric(xx))
      stop("Internal error in [make.communal]: xx not numeric")
    
    ## First, find the size of the new data frame (nn.out):
    ind.date <- cbind(xx[, 1, drop = FALSE] + xx[, 4, drop = FALSE],
                      xx[, 2, drop = FALSE] + xx[, 4, drop = FALSE])
    
    cases <- ( (ind.date[, 1] < end.per) & (ind.date[, 2] > beg.per) )
    xx <- xx[cases, , drop = FALSE]
    ind.date <- ind.date[cases, , drop = FALSE]
    
    ##  if (nrow(com.info) != n.com) stop("Error in com.info: wrong noof rows")
    ##return(xx)
    ind.iv <- get.iv(ind.date)
    ##return(ind.date, ind.iv)
    nn <- nrow(xx) ## NOTE: New definition of nn!
    nn.out <- sum(ind.iv[, 2] - ind.iv[, 1] + 1)
    ## return(nn.out)
    
    yy <- matrix(0, nrow = nn.out, ncol = ncol(xx) + 1)
    
    ## And so we fill it!
    
    nn.out <- ind.iv[, 2] - ind.iv[, 1] + 1
    cur.row <- 0

    com.dat <- as.matrix(com.dat)
    
    split <- function(i)
      {
        n.rows <- nn.out[i]
        if (n.rows == 1){
          return( c(xx[i, ], ind.iv[i, 1]) )
        }else{
          x.i <- xx[i, ]
          out <- matrix(0, nrow = n.rows, ncol = ncol(yy))
          out[n.rows, 3] <- x.i[3]
          out[, 4] <- x.i[4]
          out[, 5] <- x.i[5]
          out[1, 1] <- x.i[1]
          out[1, 2] <- cuts[ind.iv[i, 1] + 1] - x.i[4]
          out[1, 6] <- ind.iv[i, 1]
          out[n.rows, 1] <- cuts[ind.iv[i, 2]] - x.i[4]
          out[n.rows, 2] <- x.i[2]
          out[n.rows, 6] <- ind.iv[i, 2]
          if (n.rows > 2){
            for (j in 2:(n.rows - 1)){
              out[j, 1] <- out[j - 1, 2]
              out[j, 2] <- out[j - 1, 2] + iv.length
              out[j, 6] <- ind.iv[i, 1] + j - 1
            }
          }
          return(out)
        }
      }

    beg.row <- end.row <- 0

    if (!fortran){
      for (j in 1:nn){
        beg.row <- end.row + 1
        end.row <- end.row + nn.out[j]
        yy[beg.row:end.row, ] <- split(j) 
        if (j %/% 100 * 100 == j) cat("j = ", j, "\n")
        NULL
      }
    }
    if (fortran){
      
      yy <- .Fortran("split",
                     as.double(xx),
                     as.integer(nn),
                     as.integer(ncol(xx)),
                     yy = as.double(yy),
                     as.integer(nrow(yy)),
                     as.integer(ncol(yy)),
                     as.integer(nn.out),
                     as.integer(ind.iv),
                     as.double(cuts),
                     as.integer(n.years),
                     DUP = FALSE
                     )$yy
      yy <- matrix(yy, ncol = ncol(xx) + 1)
    }

    yy <- cbind(yy[, 1:4, drop = FALSE],
                dat[yy[, 5], -surv.indices, drop = FALSE],
                com.dat[yy[, 6], , drop = FALSE])
    names(yy)[1:4] <- surv
    all.names <- c(surv, fixed.names[-surv.indices], com.names)
    row.names(yy) <- as.character(1:nrow(yy))
    yy <- as.data.frame(yy)
    names(yy) <- all.names
    if (event.boolean) yy[, 3] <- as.logical(yy[, 3]) ## Beware of '3'!
  }
  
  yy
}
          
