\name{join.spells}
\alias{join.spells}

\title{Straighten up a survival data frame}
\description{
 Unnecessary cut spells are glued together, overlapping spells are
 "polished", etc. 
}
\usage{
join.spells(dat, eps = 1.e-8)
}

\arguments{
  \item{dat}{A data frame with names enter, exit, event, id.}
  \item{eps}{Tolerance for equality of two event times.}
}
\details{
 In case of overlapping intervals (i.e., a data error), the appropriate
 id's are returned.
}
\value{
  A data frame with the same variables as the input, but individual
  spells are joined, if possible (identical covariate values, and
  adjacent time intervals).
}
\references{Therneau, T.M. and Grambsch, P.M. (2000). \emph{Modeling
    Survival Data: Extending the Cox
  model.} Springer.}
\author{Göran Broström}
%\note{}


\seealso{\code{\link{coxreg}}, \code{\link{mlreg}}, \code{\link{check.surv}}}

%\examples{
%}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{survival}
