\name{cro}
\alias{cro}
\title{Creates a minimal representation of a data frame.}
\description{
Given a data frame with a defined response variable, this function
creates a unique representation of the covariates in the data frame,
vector (matrix) of responses, and a pointer vector, connecting the
responses with the corresponding covariates. 
}
\usage{
cro(dat, response=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame}
  \item{response}{The column(s) where the response resides.}
}
\details{
  The rows in the data frame are converted to text strings with
  \code{paste} and compared with \code{match}.
  
}
\value{
A list with components
  \item{y}{The response.}
  \item{covar}{A data frame with unique rows of covariates.}
  \item{keys}{Pointers from \code{y} to \code{covar}, connecting each
    response with its covariate vector.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
\note{This function is based on suggestions by Anne York and Brian Ripley.}

\seealso{\code{\link{match}}, \code{\link{paste}}}

\examples{
dat <- data.frame(y = c(1.1, 2.3, 0.7), x1 = c(1, 0, 1), x2 = c(0, 1, 0))
cro(dat)
}
\keyword{manip}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
