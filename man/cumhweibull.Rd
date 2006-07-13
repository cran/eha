\name{Hweibull}
\alias{Hweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cumulative hazard function for the Weibull distribution}
\description{
  Calculates the cumulative hazard function of a Weibull distribution 
}
\usage{
Hweibull(x, shape, scale = 1, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of quantiles.}
  \item{shape}{shape parameter.}
  \item{scale}{scale parameter, defaulting to 1.}
  \item{log}{logical: if TRUE, the answer is given as log(H).}
}
\details{
  See \link{pweibull}
}
\value{
  Gives the cumulative hazard function at x.

}
\author{Göran broström}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, shape, scale = 1, log = FALSE){
    res <- (x / scale)^shape
    if (log) res <- logb(res, base = exp(1))
    res
  }
}
\keyword{ survival }% at least one, from doc/KEYWORDS

