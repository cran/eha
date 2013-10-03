\name{Lognormal}
\alias{Lnorm}
\alias{hlnorm}
\alias{Hlnorm}
\title{The Lognormal Distribution}
\description{
  Hazard function and
  cumulative hazard function for the Lognormal distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
hlnorm(x, meanlog = 0, sdlog = 1,
shape = 1 / sdlog, scale = exp(meanlog), prop = 1, log = FALSE)
Hlnorm(x, meanlog = 0, sdlog = 1,
shape = 1 / sdlog, scale = exp(meanlog), prop = 1, log.p = FALSE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{meanlog, sdlog}{Mean and standard deviation in the distribution
    of the logarithm of a lognormal random variable.}
  \item{shape, scale}{shape and scale parameters, both defaulting to 1.}
  \item{prop}{A 'proportional hazards' parameter, see Details.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  \code{hlnorm} gives the hazard function, and
  \code{Hlnorm} gives the cumulative hazard function.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  This is a complement to the Lognormal distribution, see
  \code{\link[stats]{Lognormal}}  for further details.

  From eha version 2.3-0, the lognormal distribution is extended to a
  three-parameter family of distributions. The third parameter, 'prop',
  is a 'proportional hazards' parameter, simply multiplying the hazard
  function by its value.  For now only implemented in the hazard and
  cumulative hazards functions, but it will be introduced fully soon.
  Probably as a class of distributions with a new name.
}
\seealso{\code{\link{phreg}}, \code{\link{aftreg}}.}
\keyword{distribution}
