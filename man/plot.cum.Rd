\name{plot.cum}
\alias{plot.cum}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plots of Cumulative Hazards Functions.}
\description{
Nelson-Aalen plots.
}
\usage{
plot.cum(x, group=rep(1, length(exit)), 
main="Survivor function(s)", xlab="Duration", ylab="Remaining fraction",
log.scale = FALSE, ...)
}
\arguments{
  \item{x}{A \code{Surv} object.}
  \item{group}{Defines a partition of the data. One survivor function
    for each level of \code{group} is drawn.}
  \item{main}{Header of the plot.}
  \item{xlab}{Label of the x axis.}
  \item{ylab}{Label of the y axis.}
  \item{log.scale}{If TRUE, plotting is done on a log-log scale. Default
    is FALSE}
  \item{...}{Anything that \code{plot} likes...}
}
\details{
Left truncation is allowed. Note, though, that this fact may result in
strange estimated curves due to lack of data in certain (low) ages.
}
\value{
No value is returned.
}
\references{}
\author{Göran Broström}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{plot.Surv}}}

\examples{
time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50)
plot.cum(Surv(time0, time1, event), group)
plot.cum(Surv(time0, time1, event), group, log.scale = TRUE)
}
\keyword{survival}


