age.window <- function(dat, window,
                       surv = c("enter", "exit", "event")){

  if (!is.data.frame(dat))stop("dat must be a data frame")
  if (length(surv) != 3) stop("surv must have length 3")
  fixed.names <- names(dat)
  surv.indices <- match(surv, fixed.names)
  if (length(which(is.na(surv.indices)))){
    x <- which(is.na(surv.indices))
    stop(paste(surv[x], " is not a name in the data frame."))
  }

  enter <- dat[, surv.indices[1]]
  exit <- dat[, surv.indices[2]]
  event <- dat[, surv.indices[3]]
  
  event <- ifelse( (exit > window[2]), 0, event)
  who <- (exit > window[1]) & (enter < window[2])
  event <- event[who]
  exit <- pmin(exit[who], window[2])
  enter <- pmax(enter[who], window[1])

  dat <- dat[who, ]
  dat[, surv.indices[1]] <- enter
  dat[, surv.indices[2]] <- exit
  dat[, surv.indices[3]] <- event
  
  dat
}
