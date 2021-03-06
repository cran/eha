% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perstat.R
\name{perstat}
\alias{perstat}
\title{Period statistics}
\usage{
perstat(surv, period, age = c(0, 200))
}
\arguments{
\item{surv}{An (extended) \code{surv} object (4 columns with \code{enter},
\code{exit}, \code{event}, \code{birthdate})}

\item{period}{A vector of dates (in decimal form)}

\item{age}{A vector of length 2; lowest and highest age}
}
\value{
A list with components \item{events}{No. of events in eavh time
period.} \item{exposure}{Exposure times in each period.}
\item{intensity}{\code{events / exposure}}
}
\description{
Calculates occurrence / exposure rates for time periods given by
\code{period} and for ages given by \code{age}.
}
\seealso{
\code{\link{piecewise}}
}
\author{
Göran Broström
}
\keyword{nonparametric}
\keyword{survival}
