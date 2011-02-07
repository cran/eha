\name{join.spells}
\alias{join.spells}

\title{Straighten up a survival data frame}
\description{
 Unnecessary cut spells are glued together, overlapping spells are
 "polished", etc. 
}
\usage{
join.spells(dat, strict = FALSE, eps = 1.e-8)
}

\arguments{
  \item{dat}{A data frame with names enter, exit, event, id.}
  \item{strict}{If TRUE, nothing is changed if errors in spells
    (non-positive length, overlapping intervals, etc.) are
    detected. Otherwise (the default), bad spells are removed, with
    "earlier life" having higher priority.}
  \item{eps}{Tolerance for equality of two event times. Should be kept
    small.}
}
\details{
 In case of overlapping intervals (i.e., a data error), the appropriate
 id's are returned if \code{strict} is \code{TRUE}.
}
\value{
  A data frame with the same variables as the input, but individual
  spells are joined, if possible (identical covariate values, and
  adjacent time intervals).
}
\references{Therneau, T.M. and Grambsch, P.M. (2000). \emph{Modeling
    Survival Data: Extending the Cox
  model.} Springer.}
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{}


\seealso{\code{\link{coxreg}}, \code{\link{aftreg}}, \code{\link{check.surv}}}

%\examples{
%}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{survival}
