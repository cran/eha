\name{plot.hazdata}
\alias{plot.hazdata}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plots of survivor functions.}
\description{
  Baseline hazards estimates.
}
\usage{
plot.hazdata(x, fn = c("cum", "surv", "log", "loglog"), fig = TRUE,
xlim=NULL, ylim=NULL, main=NULL, xlab=NULL, ylab=NULL, ...)
}
\arguments{
  \item{x}{A \code{hazdata} object, typically the 'hazards' element in
the output from \code{link{coxreg}}.}
  \item{fn}{Which type of plot?}
  \item{fig}{Should a plot actually be produced? Default is TRUE.}
  \item{xlim}{Horizontal plot limits. If NULL, calculated by the function.}
  \item{ylim}{Vertical plot limits. If NULL, set to \code{c(0, 1)}}
  \item{main}{A heading for the plot.}
  \item{xlab}{Label on the x axis.}
  \item{ylab}{Label on the y-axis.}
  \item{...}{Anything that \code{\link[base]{plot.default}} likes...}
}
\details{
It is also possible to have as first argument an object of type
"coxreg", given that it contains a component of type "hazdata".
}
\value{
A list where the elements are two-column matrices, one for each stratum
in the model. The first column contains risktimes, and the second the y
coordinates for the requested curve(s).
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
plot(fit$hazards)
}
\keyword{survival}
