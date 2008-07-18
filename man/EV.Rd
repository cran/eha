\name{EV}
\alias{EV}
\alias{dEV}
\alias{pEV}
\alias{qEV}
\alias{hEV}
\alias{HEV}
\alias{rEV}
\title{The EV Distribution}
\description{
  Density, distribution function, quantile function, hazard function,
  cumulative hazard function, and random
  generation for the EV distribution with parameters \code{shape}
  and \code{scale}.
}
\usage{
dEV(x, shape = 1, scale = 1, log = FALSE)
pEV(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
qEV(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
hEV(x, shape = 1, scale = 1, log = FALSE)
HEV(x, shape = 1, scale = 1, log.p = FALSE)
rEV(n, shape = 1, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape, scale}{shape and scale parameters, both defaulting to 1.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dEV} gives the density,
  \code{pEV} gives the distribution function,
  \code{qEV} gives the quantile function,
  \code{hEV} gives the hazard function,
  \code{HEV} gives the cumulative hazard function, and
  \code{rEV} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  The EV distribution with \code{shape} parameter \eqn{a} and
  \code{scale} parameter \eqn{\sigma}{b} has hazard given by
  \deqn{h(x) = \exp((x / \sigma)^b)}{%
        f(x) = exp((x / a)^b)}
  for \eqn{x \ge 0}{x >= 0}.
}
\keyword{distribution}
