\name{hweibull}
\alias{hweibull}

\title{The hazard function of a Weibull distribution}
\description{
 Calculates the hazard function of a Weibull distribution. 
}
\usage{
hweibull(x, shape, scale = 1, log = FALSE)
}

\arguments{
  \item{x}{vector of quantiles}
  \item{shape}{The shape parameter}
  \item{scale}{The scale parameter, defaults to 1.}
  \item{log}{logical; if TRUE, the log of the hazard function is given.}
}
\details{
  See \link{dweibull}.
}
\value{
The hazard function, evaluated at x.
}
\author{Göran Broström}
\seealso{\code{\link{pweibull}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, shape, scale = 1, log = FALSE){
    if (shape <= 0 || scale <= 0)
      error("scale and shape must be positive")  
    (x / scale)^(shape - 1) / scale
  }
}
\keyword{ survival }% at least one, from doc/KEYWORDS
