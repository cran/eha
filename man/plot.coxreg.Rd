\name{plot.coxreg}
\alias{plot.coxreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plots of survivor functions.}
\description{
  Baseline hazards estimates.
}
\usage{
\method{plot}{coxreg}(x, fn = c("cum", "surv", "log", "loglog"), fig = TRUE,
xlim=NULL, ylim=NULL, main=NULL, xlab="Duration", ylab="",
new.data = NULL, \ldots)
}
\arguments{
  \item{x}{A \code{coxreg} object, typically the output from
    \code{link{coxreg}}.} 
  \item{fn}{Which type of plot?}
  \item{fig}{Should a plot actually be produced? Default is TRUE.}
  \item{xlim}{Horizontal plot limits. If NULL, calculated by the function.}
  \item{ylim}{Vertical plot limits. If NULL, set to \code{c(0, 1)} for a
  survival plot, otherwise adaptive, data dependent.}
  \item{main}{A heading for the plot.}
  \item{xlab}{Label on the x axis.}
  \item{ylab}{Label on the y-axis.}
  \item{new.data}{At what covariate values should the calculations be
    done? Default is the mean values of the covariates.}
  \item{...}{Anything that \code{\link[graphics]{plot.default}} likes...}
}
\details{
This function is a wrapper for \code{\link{plot.hazdata}}.
}
\value{
A list where the elements are two-column matrices, one for each stratum
in the model. The first column contains risktimes, and the second the y
coordinates for the requested curve.
}
%\references{}
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\examples{
time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50) 
fit <- coxreg(Surv(time0, time1, event) ~ strata(group))
plot.coxreg(fit)
}
\keyword{survival}
