\name{table.events}
\alias{table.events}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculating failure times, risk set sizes and No. of events in
  each risk set}
\description{
  From input data of the 'interval' type, with an event indicator,
  summary statistics for each risk set (at an event time point) are
  calculated. 
}
\usage{
table.events(enter=rep(0, length(exit)), exit, event, strict=TRUE )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{enter}{Left truncation time point.}
  \item{exit}{End time point, an event or a right censoring.}
  \item{event}{Event indicator.}
  \item{strict}{If TRUE, then tabulating is not done after a time point
    where all individuals in a riskset failed.}
}
%\details{
  
%}
\value{
  A list with components
  \item{times}{Ordered distinct event time points.}
  \item{events}{Number of events at each event time point.}
  \item{riskset.sizes}{Number at risk at each event time point.}
}
%\references{}
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{}

 
\seealso{\code{\link{risksets}}}

\examples{
exit = c(1,2,3,4,5)
event = c(1,1,0,1,1)
table.events(exit = exit, event = event)
}
\keyword{survival}%-- one or more ...
