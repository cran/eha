cro <- function(dat, response = 1){
  covar <- unique(dat[, -response, drop = FALSE])
  dat.keys <-
    match(do.call("paste", c(dat[, -response, drop = FALSE], sep="\001")),
          do.call("paste", c(covar,  sep="\001")))

  return(y = dat[, response],
         covar = covar,
         keys = dat.keys)
}
