\name{Gompertz}
\alias{Gompertz}
\alias{dgompertz}
\alias{pgompertz}
\alias{qgompertz}
\alias{hgompertz}
\alias{Hgompertz}
\alias{rgompertz}
\title{The Gompertz-Gompertz Distribution}
\description{
  Density, distribution function, quantile function, hazard function,
  cumulative hazard function, and random
  generation for the Gompertz-Gompertz distribution with parameters
  \code{shape1}, \code{shape2}
  and \code{scale}.
}
\usage{
dgompertz(x, shape = 1, scale = 1, log = FALSE)
pgompertz(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
qgompertz(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
hgompertz(x, shape = 1, scale = 1, log = FALSE)
Hgompertz(x, shape = 1, scale = 1, log.p = FALSE)
rgompertz(n, shape = 1, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape, scale}{Parameters: shape, defaulting to 1, 
    and scale, defaulting to 1.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dgompertz} gives the density,
  \code{pgompertz} gives the distribution function,
  \code{qgompertz} gives the quantile function,
  \code{hgompertz} gives the hazard function,
  \code{Hgompertz} gives the cumulative hazard function, and
  \code{rgompertz} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  The Gompertz distribution with \code{shape} parameter \eqn{a} and
  \code{scale} parameter \eqn{\sigma}{b} has hazard given by
  \deqn{h(x) = a \exp(x/\sigma)}{%
        h(x) = a exp(x/b)}
  for \eqn{x \ge 0}{x >= 0}.
}
\keyword{distribution}
