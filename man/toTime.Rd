\name{toTime}
\alias{toTime}

\title{Calculate duration in years from "0000-01-01" to a given date}
\description{
 Given a vector of dates, the output is a vector of durations in years
 since "0000-01-01".
}
\usage{
toTime(dates)
}
\arguments{
  \item{dates}{A vector of dates in character form or of class \code{Date}}
}

\value{
A vector of durations, as decribed above.
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{ ~~further notes~~ }

\seealso{\code{\link{toDate}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
toTime(c("1897-05-16", "1901-11-21"))
}

\keyword{survival}% at least one, from doc/KEYWORDS

