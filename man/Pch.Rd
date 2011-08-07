\name{pch}
\alias{pch}
\alias{dpch}
\alias{ppch}
\alias{qpch}
\alias{hpch}
\alias{Hpch}
\alias{rpch}
\title{The Piecewise constant hazards (Pch) distribution}
\description{
  Density, distribution function, quantile function, hazard function,
  cumulative hazard function, and random
  generation for the pch distribution with parameters \code{cuts}
  and \code{levels}.
}
\usage{
dpch(x, cuts, levels, log = FALSE)
ppch(q, cuts, levels, lower.tail = TRUE, log.p = FALSE)
qpch(p, cuts, levels, lower.tail = TRUE, log.p = FALSE)
hpch(x, cuts, levels, log = FALSE)
Hpch(x, cuts, levels, log.p = FALSE)
rpch(n, cuts, levels)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{cuts, levels}{\code{cuts} devine the intervals where the hazard
    function is constant. The cuts must be strictly positive and
    finite. \code{levels} are the interval-constant values, one more
    than the cuts.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.}
}
\value{
  \code{dpch} gives the density,
  \code{ppch} gives the distribution function,
  \code{qpch} gives the quantile function,
  \code{hpch} gives the hazard function,
  \code{Hpch} gives the cumulative hazard function, and
  \code{rpch} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\details{
  The Pch distribution is defined by the cuts and the levels so that the
  hazard function is constant on intervals.
  
}
\keyword{distribution}
