%\encoding{UTF-8}
\name{check.surv}
\alias{check.surv}
\title{Check the integrity of survival data.}
\description{
Check that exit occurs after enter, that spells from an individual do
not overlap, and that each individual experiences at most one event.
}
\usage{
check.surv(enter, exit, event, id = NULL, eps = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enter}{Left truncation time.}
  \item{exit}{Time of exit.}
  \item{event}{Indicator of event. Zero means 'no event'.}
  \item{id}{Identification of individuals.}
  \item{eps}{The smallest allowed spell length or overlap.}
}
\details{
Interval lengths must be strictly positive.
}
\value{
A vector of id's for the insane individuals. Of zero length if no errors.
}
\references{}
\author{Göran Broström}
\note{}

\seealso{\code{\link{join.spells}}, \code{\link{coxreg}}, \code{\link{mlreg}}}
\examples{
xx <- data.frame(enter = c(0, 1), exit = c(1.5, 3), event = c(0, 1), id =
c(1,1))
check.surv(xx$enter, xx$exit, xx$event, xx$id)
}
\keyword{manip}
\keyword{survival}
