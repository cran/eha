\name{clean}
\alias{clean}

\title{Clean a survival data frame}
\description{
 Unnecessary cut spells are glued together, overlapping spells are
 "polished", etc. 
}
\usage{
clean(dat)
}

\arguments{
  \item{dat}{A data frame with names enter, exit, event, id, birthdate. }
}
\details{
 In case of overlapping intervals (i.e., a data error), "reasonable"
 changes are made. Everyone may not agree upon what is reasonable,
 though, so check the code to see exactly what happens.  
}
\value{
  A data frame with the same variables as the input.
}
\references{Therneau, T.M. and Grambsch, P.M. (2000). \emph{Modeling
    Survival Data: Extending the Cox
  model.} Springer.}
\author{Göran Broström}
%\note{}


\seealso{\code{\link{coxreg}}, \code{\link{mlreg}}, }

%\examples{
%}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{survival}
