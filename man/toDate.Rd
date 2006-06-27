%\encoding{UTF-8}
\name{toDate}
\alias{toDate}

\title{Convert time in years since "0000-01-01" to a date.}
\description{
  This function uses \code{as.Date} and a simple linear transformation.
}
\usage{
toDate(times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{a vector of durations}
}
\details{

}
\value{
A vector of dates as character strings of the type "1897-05-21".
}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
%\note{ ~~further notes~~ }


\seealso{ \code{\link{toTime}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
toDate(1897.357)
}
\keyword{survival}% at least one, from doc/KEYWORDS
