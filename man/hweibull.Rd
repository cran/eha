\encoding{UTF-8}
\name{hweibull}
\alias{hweibull}
\alias{Hweibull}

\title{The (Cumulative) Hazard Function of a Weibull Distribution}
\description{
 \code{hweibull} calculates the hazard function of a Weibull
 distribution, and \code{Hweibull} calculates the corresponding
 cumulative hazard function. 
}
\usage{
hweibull(x, shape, scale = 1, log = FALSE)
Hweibull(x, shape, scale = 1, log = FALSE)
}

\arguments{
  \item{x}{Vector of quantiles.}
  \item{shape}{The shape parameter.}
  \item{scale}{The scale parameter, defaults to 1.}
  \item{log}{logical; if TRUE, the log of the hazard function is given.}
}
\details{
  See \link{dweibull}.
}
\value{
The (cumulative) hazard function, evaluated at x.
}
\author{\enc{Göran Broström}{George Bridgewater}}
\seealso{\code{\link{pweibull}}}
\examples{
hweibull(3, 2, 1)
dweibull(3, 2, 1) / pweibull(3, 2, 1, lower.tail = FALSE)
Hweibull(3, 2, 1)
-pweibull(3, 2, 1, log.p = TRUE, lower.tail = FALSE)
}
\keyword{ survival }% at least one, from doc/KEYWORDS
