%\encoding{UTF-8}
\name{frail.fit}
\alias{frail.fit}

\title{Fits a frailty proportional hazards model}
\description{
This function is called from \code{\link{mlreg}}. Can be called directly
by a user who knows what (s)he is doing! 
}
\usage{
frail.fit(X, Y, rs, strats, offset, init, max.survs, frailty, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Design matrix}
  \item{Y}{Survival object}
  \item{rs}{Risk set pointers}
  \item{strats}{Stratum variable}
  \item{offset}{Offset variable}
  \item{init}{Initial regression parameter values}
  \item{max.survs}{Sampling of risk sets?}
  \item{frailty}{The frailty grouping variable}
  \item{control}{Values for the optimization}
}
\details{
Frailty is analyzed in a marginal framework, with numerical
integration. This is very slow and somewhat unstable.
}
\value{
See the code!
}

\author{Göran Broström}

\section{Warning }{Not very much tested. Can fail, and is very slow.}
\seealso{\code{\link{mlreg}}}
\examples{
fam.size <- rep(20, 5)
n <- sum(fam.size)
fam <- rep(1:length(fam.size), fam.size)
beta <- 1
sigma <- 1
frail <- rep(rnorm(length(fam.size), 0, sigma), fam.size)
x <- runif(n)
exit <- rexp(n, exp(x * beta + frail))
event <- rep(1, n)
enter <- numeric(n)
mlreg(Surv(enter, exit, event) ~ x, frailty = frail)
}
\keyword{survival}
\keyword{regression}
