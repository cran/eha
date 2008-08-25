\name{Makeham}
\alias{Makeham}
\alias{dmakeham}
\alias{pmakeham}
\alias{qmakeham}
\alias{hmakeham}
\alias{Hmakeham}
\alias{rmakeham}
\title{The Gompertz-Makeham Distribution}
\description{
  Density, distribution function, quantile function, hazard function,
  cumulative hazard function, and random
  generation for the Gompertz-Makeham distribution with parameters
  \code{shape1}, \code{shape2}
  and \code{scale}.
}
\usage{
dmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
pmakeham(q, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
qmakeham(p, shape = c(1, 1), scale = 1, lower.tail = TRUE, log.p = FALSE)
hmakeham(x, shape = c(1, 1), scale = 1, log = FALSE)
Hmakeham(x, shape = c(1, 1), scale = 1, log.p = FALSE)
rmakeham(n, shape = c(1, 1), scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape, scale}{Parameters: shape, a vector of length 2,
    defaulting to c(1, 1), 
    and scale, defaulting to 1.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}
}
\value{
  \code{dmakeham} gives the density,
  \code{pmakeham} gives the distribution function,
  \code{qmakeham} gives the quantile function, but is not yet implemented,
  \code{hmakeham} gives the hazard function,
  \code{Hmakeham} gives the cumulative hazard function, and
  \code{rmakeham} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  The Gompertz-Makeham distribution with \code{shape} parameters
  \eqn{a_1}{a1} and \eqn{a_2}{a2} and
  \code{scale} parameter \eqn{\sigma}{b} has hazard function given by
  \deqn{h(x) = a_2 + a_1 \exp(x / \sigma)}{%
        h(x) = a2 + a1 exp(x / b)}
  for \eqn{x \ge 0}.
}
\keyword{distribution}
