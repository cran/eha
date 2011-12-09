\name{plot.aftreg}
\alias{plot.aftreg}

\title{Plots output from a Weibull regression}
\description{
Just a simple plot of the hazard functions for each stratum.
}
\usage{
\method{plot}{aftreg}(x, fn = c("haz", "cum", "den", "sur"), main = NULL,
xlim = NULL, ylim = NULL, xlab = "Duration", ylab = "",
new.data = x$means, \ldots)
}

\arguments{
  \item{x}{A \code{aftreg} object}
  \item{fn}{Which functions shoud be plotted! Default is all. They
    will scroll by, so you have to take care explicitely what you want
    to be produced. See, eg, \code{par(mfrow = ...)}}
  \item{main}{Header for the plot}
  \item{xlim}{x limits}
  \item{ylim}{y limits}
  \item{xlab}{x label}
  \item{ylab}{y label}
  \item{new.data}{At which covariate values?}
  
  \item{\dots}{Extra parameters passed to 'plot'}
}

\details{
The plot is drawn at the mean values of the covariates, by default.
}
\value{
No return value.
}

\author{\enc{Göran Broström}{George Bridgewater}}

\seealso{\code{\link{aftreg}}}
\examples{
y <- rllogis(40, shape = 1, scale = 1)
x <- rep(c(1,1,2,2), 10)
fit <- aftreg(Surv(y, rep(1, 40)) ~ x, dist = "loglogistic")
plot(fit)
}
\keyword{dplot}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
