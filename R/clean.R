## A very special function for the 'Early Life' project. Göran Broström

clean <- function(dat){
  eps <- 0.0000001
  ## Assumes: enter, exit, event, id, birthdate
  #Must have original id (as a covariate):
  resp <- match(c("enter", "exit", "event",
                  "startevent", "stopevent" ), names(dat))
  if (any(is.na(resp))) stop("Wrong variable names")

  covar <- dat[ , -resp]
  n.cov <- ncol(covar)
  n.rec <- nrow(covar)
  all <- unique(dat$id)
  nn <- length(all)

  ide <- as.integer(factor(dat$id, labels = 1:nn))
  ord <- order(ide, dat$enter, dat$exit)
  dat <- dat[ord, ]

  res <- .Fortran("cleanup",
                  as.double(t(covar)),
                  as.double(dat$enter),
                  as.double(dat$exit),
                  as.integer(dat$event),
                  as.integer(ide),
                  as.integer(n.cov),
                  as.integer(n.rec),
                  as.integer(nn),
                  new.n.rec = integer(1),
                  new.cov = double(n.rec * n.cov),
                  enter = double(n.rec),
                  exit = double(n.rec),
                  event = integer(n.rec),
                  id = integer(n.rec),
                  ##DUP = FALSE,
                  PACKAGE = "eha")

  
  out <- data.frame(new.id = res$id[1:res$new.n.rec],
                    enter = res$enter[1:res$new.n.rec],
                    exit = res$exit[1:res$new.n.rec],
                    event = res$event[1:res$new.n.rec]
                    )

  new.cov <- data.frame(matrix(res$new.cov, byrow = TRUE,
                               ncol = n.cov))[1:res$new.n.rec, , drop = FALSE]

  names(new.cov) <- names(covar)

  cbind(out, new.cov)
}
