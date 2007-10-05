\name{perstat}
\alias{perstat}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Period statistics}
\description{
  Calculates occurrence / exposure rates for time periods given by
  \code{period} and for ages given by \code{age}.
  
}
\usage{
perstat(surv, period, age = c(0, 200))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{surv}{An (extended) \code{surv} object (4 columns with
    \code{enter}, \code{exit}, \code{event}, \code{birthdate}) }
  \item{period}{A vector of dates (in decimal form)}
  \item{age}{A vector of length 2; lowest and highest age}
}
\details{
 
}
\value{
A list with components
  \item{events}{No. of events in eavh time period.}
  \item{exposure}{Exposure times in each period.}
  \item{intensity}{\code{events / exposure}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{ ~~further notes~~ }

\seealso{\code{\link{piecewise}}}

%\examples{
%}
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{nonparametric}% __ONLY ONE__ keyword per line
