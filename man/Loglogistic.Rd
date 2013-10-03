\name{Loglogistic}
\alias{Llogis}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{hllogis}
\alias{Hllogis}
\alias{rllogis}
\title{The Loglogistic Distribution}
\description{
  Density, distribution function, quantile function, hazard function,
  cumulative hazard function, and random
  generation for the Loglogistic distribution with parameters \code{shape}
  and \code{scale}.
}
\usage{
dllogis(x, shape = 1, scale = 1, log = FALSE)
pllogis(q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
qllogis(p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
hllogis(x, shape = 1, scale = 1, prop = 1, log = FALSE)
Hllogis(x, shape = 1, scale = 1, prop = 1, log.p = FALSE)
rllogis(n, shape = 1, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape, scale}{shape and scale parameters, both defaulting to 1.}
  \item{prop}{A 'proportional hazards' parameter, for now only available
  for \code{hllogis} and \code{Hllogis}. See Details.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dllogis} gives the density,
  \code{pllogis} gives the distribution function,
  \code{qllogis} gives the quantile function,
  \code{hllogis} gives the hazard function,
  \code{Hllogis} gives the cumulative hazard function, and
  \code{rllogis} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  The Loglogistic distribution with \code{shape} parameter \eqn{a} and
  \code{scale} parameter \eqn{\sigma}{b} has density given by
  \deqn{f(x) = (a/\sigma) {(x/\sigma)}^{a-1} / (1 + {(x/\sigma)}^{a})^2}{%
        f(x) = (a/b) (x/b)^(a-1) / (1 + (x/b)^a)^2}
  for \eqn{x \ge 0}.
  The cumulative distribution function is
  \eqn{F(x) = 1 - 1 / (1 + {(x/\sigma)}^a)}{F(x) = 1 - 1 / (1 + (x/b)^a)}
  on \eqn{x \ge 0}.

  From eha version 2.3-0, the loglogistic distribution is extended to a
  three-parameter family of distributions. The third parameter, 'prop',
  is a 'proportional hazards' parameter, simply multiplying the hazard
  function by its value.  For now only implemented in the hazard and
  cumulative hazards functions, but it will be introduced fully soon.
  Probably as a class of distributions with a new name.

  See also \code{\link{phreg}}.
}
\seealso{\code{\link{phreg}}, \code{\link{aftreg}}.}
\keyword{distribution}
