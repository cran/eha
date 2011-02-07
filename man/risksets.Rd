\name{risksets}
\alias{risksets}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Finds the compositions and sizes of risk sets}
\description{
  Focus is on the risk set composition just prior to a failure.
}
\usage{
risksets(x, strata, max.survs)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A \code{Surv} object.}
  \item{strata}{Stratum indicator.}
  \item{max.survs}{Maximum number of survivors in each risk set. If
    smaller than the 'natural number', survivors are sampled from the
    present ones.}
}
\details{
  If the input argument max.survs is left alone, all survivors are
  accounted for in all risk sets.
}
\value{
  A list with components
  \item{antrs}{No. of risk sets in each stratum. The number of strata is
    given by \code{length(antrs)}.}
  \item{risktimes}{Ordered distinct failure time points.}
  \item{eventset}{vector of pointers to events in each risk set.} 
  \item{riskset}{vector of pointers to 
    the members of the risk sets, in order. The 'n.events' 
    first are the events.}
  \item{size}{The sizes of the risk sets.}
  \item{n.events}{The number of events in each risk set.}
}

\author{\enc{Göran Broström}{George Bridgewater}}
\note{can be used to "sample the risk sets".}

\seealso{\code{\link{table.events}}, \code{\link{coxreg}}.}

\examples{
 enter <- c(0, 1, 0, 0)
 exit <- c(1, 2, 3, 4)
 event <- c(1, 1, 1, 0)
 risksets(Surv(enter, exit, event))
}
\keyword{survival}
