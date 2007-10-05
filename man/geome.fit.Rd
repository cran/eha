\name{geome.fit}
\alias{geome.fit}

\title{Constant intensity discrete time proportional hazards}
\description{
This function is called from \code{mlreg}. A user may call it directly.
}
\usage{
geome.fit(X, Y, rs, strats, offset, init, max.survs, method = "ML",
boot = FALSE, control)
}

\arguments{
  \item{X}{The design matrix}
  \item{Y}{Survival object}
  \item{rs}{risk set produced by \code{risksets}}
  \item{strats}{Stratum indicator}
  \item{offset}{ Offset}
  \item{init}{Initial values}
  \item{max.survs}{Maximal survivors}
  \item{method}{"ML" or "MPPL"}
  \item{boot}{should we bootstrap?}
  \item{control}{See mlreg}
}

\value{
See the code.
}
\references{See mlreg}
\author{Göran Broström}
\note{Nothing special}

\seealso{\code{\link{mlreg}}}
\examples{

}
\keyword{survival}

