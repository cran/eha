\name{plot.Surv}
\alias{plot.Surv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plots of survivor functions.}
\description{
  Kaplan-Meier estimates. If only one curve, confidence limits according
  to Greenwood's formula are drawn.
}
\usage{
plot.Surv(x, strata=rep(1, length(exit)), limits=FALSE, conf=0.95,
xlim = NULL, ylim = NULL, main="Survivor function(s)", 
xlab="Duration",
ylab="Remaining fraction", ...)
}
\arguments{
  \item{x}{A \code{Surv} object.}
  \item{strata}{Defines a partition of the data. One survivor function
    for each level of \code{strata} is drawn.}
  \item{limits}{If TRUE, and if the number of curves is one, confidence
    limits are drawn.}
  \item{conf}{The confidence level for the confidence limits.}
  \item{xlim}{Horizontal plot limits. If NULL, calculated by the function.}
  \item{ylim}{Vertical plot limits. If NULL, set to \code{c(0, 1)}}.
  \item{main}{A heading for the plot.}
  \item{xlab}{Label on the x axis.}
  \item{ylab}{Label on the y-axis.}
  \item{...}{Anything that \code{plot} likes...}
}
\details{
Left truncation is allowed. Note, though, that this fact may result in
strange estimated curves due to lack of data in certain (low) ages.
}
\value{
No value is returned.
}
%\references{}
\author{Göran Broström}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{plot.cum}}}

\examples{
time0 <- numeric(50)
group <- c(rep(0, 25), rep(1, 25))
time1 <- rexp( 50, exp(group) )
event <- rep(1, 50)
plot.Surv(Surv(time0, time1, event), strata = group)
}
\keyword{survival}
