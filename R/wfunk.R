# 'For internal use'

wfunk <- function(beta, lambda, p, X, Y,
                  offset = rep(0, NROW(X)),
                  ord = 2, pfixed = FALSE,
                    trace = TRUE){

## Returns loglik, score, and information (=-fpp) 
## For one stratum (only)!!

  nn <- NROW(Y)
  if (NCOL(Y) == 2) Y <- cbind(rep(0, nn), Y)
  mb <- NCOL(X)
  if (pfixed){
    bdim <- mb + 1
    b <- c(beta, log(lambda))
  }else{
    bdim <- mb + 2
    b <- c(beta, log(lambda), log(p))
  }
    
  fit <- .Fortran("wfunc", ## Returns -loglik, -score, +information
                  as.integer(ord),
                  as.integer(pfixed),
                  as.double(p),
                  as.integer(bdim),
                  as.integer(mb),
                  as.double(b),
                  #
                  as.integer(nn),
                  as.double(t(X)),
                  as.double(Y[, 1]), ## enter
                  as.double(Y[, 2]), ## exit
                  as.integer(Y[, 3]), ## event
                  as.double(offset),
                  #
                  f = double(1),
                  fp = double(bdim),
                  fpp = double(bdim * bdim), ok = integer(1),
                  PACKAGE = "eha")

  if (ord >= 2) fit$fpp <- matrix(fit$fpp, ncol = bdim)
  if (ord <= 1) fit$fpp <- NULL
  if (ord <= 0) fit$fp <- NULL
  list(f = -fit$f,
       fp = -fit$fp,
       fpp = fit$fpp)
}
           
