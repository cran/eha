\name{piecewise}
\alias{piecewise}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise hazards}
\description{
Calculate piecewise hazards, no. of events, and exposure times in each
interval indicated by cutpoints.
}
\usage{
piecewise(enter, exit, event, cutpoints)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{enter}{Left interval endpoint}
  \item{exit}{Right interval endpoint}
  \item{event}{Indicator of event}
  \item{cutpoints}{Vector of cutpoints}
}
\details{
Exact calculation.
}
\value{
  A list with components

  \item{deaths}{Vector of number of deaths}
  \item{exposure}{Vector of total exposure time}
  \item{hazard}{Vector of hazards, \code{hazard} == \code{deaths} /
    \code{exposure}} 
}

\author{Göran Broström}

\seealso{\code{\link{perstat}}}

%\examples{
%}
\keyword{survival}
\keyword{nonparametric}

