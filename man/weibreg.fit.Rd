\name{weibreg.fit}
\alias{weibreg.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weibull regression}
\description{
This function is called by \code{\link{weibreg}}, but it can also be
directly called by a user.
}
\usage{
weibreg.fit(X, Y, strata, offset, init, shape, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The design (covariate) matrix.}
  \item{Y}{A survival object, the response.}
  \item{strata}{A stratum variable.}
  \item{offset}{Offset.}
  \item{init}{Initial regression parameter values.}
  \item{shape}{If positive, a fixed value of the shape parameter in the
    Weibull distribution. Otherwise, the shape is estimated.}
  \item{control}{Controls convergence and output.}
}
\details{
See \code{\link{weibreg}} for more detail.
}
\value{
  \item{coefficients}{Estimated regression coefficients plus estimated
    scale and shape coefficients, sorted by strata, if present.}
  \item{var}{}
  \item{loglik}{Vector of length 2. The first component is the maximized
  loglihood with only scale and shape in the model, the second the final
  maximum.}
  \item{score}{Score test statistic at initial values}
  \item{linear.predictors}{Linear predictors for each interval.}
  \item{means}{Means of the covariates}
  \item{conver}{TRUE if convergence}
  \item{fail}{TRUE if failure}
  \item{iter}{Number of Newton-Raphson iterates.}
  \item{n.strata}{The number of strata in the data.}

}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
%\note{}


\seealso{\code{\link{weibreg}}}
%\examples{

%}
\keyword{survival}
\keyword{regression}
