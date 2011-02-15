\name{check.dist}
\alias{check.dist}

\title{Graphical goodness-of-fit test}
\description{
Comparison of the cumulative hazards functions for a semi-parametric and
a parametric model.
}
\usage{
check.dist(sp, pp, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{An object of type "cotxreg", typically output from
\code{\link{coxreg}}} 
  \item{pp}{An object of type "phreg", typically output from
\code{\link{phreg}}}
\item{main}{Header for the plot. Default is distribution and "cumulativa
  hazard function"}
}
\details{
  For the moment only a graphical comparison.
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
fit.ln <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "lognormal")
fit.ev <- phreg(Surv(enter, exit, event) ~ ses, data = mort,
dist = "ev")
check.dist(fit.cr, fit.w)
check.dist(fit.cr, fit.g)
check.dist(fit.cr, fit.ln)
check.dist(fit.cr, fit.ev)
par(oldpar)
}
\keyword{distribution}

