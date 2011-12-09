\name{plot.weibreg}
\alias{plot.weibreg}

\title{Plots output from a Weibull regression}
\description{
Just a simple plot of the hazard functions for each stratum.
}
\usage{
\method{plot}{weibreg}(x, fn = c("haz", "cum", "den", "sur"), main = NULL,
xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
new.data = x$means, \ldots)
}

\arguments{
  \item{x}{A \code{weibreg} object}
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
The plot is drawn at the mean values of the covariates.
}
\value{
No return value
}

\author{\enc{Göran Broström}{George Bridgewater}}

\seealso{\code{\link{weibreg}}}
\examples{
y <- rweibull(4, shape = 1, scale = 1)
x <- c(1,1,2,2)
fit <- weibreg(Surv(y, c(1,1,1,1)) ~ x)
plot(fit)
}
\keyword{dplot}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
