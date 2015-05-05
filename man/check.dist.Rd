\name{check.dist}
\alias{check.dist}

\title{Graphical goodness-of-fit test}
\description{
Comparison of the cumulative hazards functions for a semi-parametric and
a parametric model.
}
\usage{
check.dist(sp, pp, main = NULL, col = NULL, printLegend = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{An object of type "coxreg", typically output from
\code{\link{coxreg}}} 
  \item{pp}{An object of type "phreg", typically output from
\code{\link{phreg}}}
\item{main}{Header for the plot. Default is distribution and "cumulative
  hazard function"}
\item{col}{Line colors. should be \code{NULL} (black lines) or of length
  2}
\item{printLegend}{Should a legend be printed? Default is \code{TRUE}.}
}
\details{
  For the moment only a graphical comparison. The arguments \code{sp}
  and \code{pp} may be swapped.
}
\value{
No return value.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}

\seealso{\code{\link{coxreg}} and \code{\link{phreg}}.}
\examples{
data(mort)
oldpar <- par(mfrow = c(2, 2))
fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
fit.w <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
fit.g <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "gompertz")
fit.pch <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "pch", cuts = c(5, 10, 15))
fit.ev <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "ev")
check.dist(fit.cr, fit.w)
check.dist(fit.cr, fit.g)
check.dist(fit.cr, fit.pch)
check.dist(fit.cr, fit.ev)
par(oldpar)
}
\keyword{distribution}

