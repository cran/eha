\name{wfunk}
\alias{wfunk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Loglihood function of a Weibull regression}
\description{
Calculates minus the log likelihood function and its first and second order
derivatives for data from a Weibull regression model. Is called by
\code{\link{weibreg}}. 
}
\usage{
wfunk(beta, lambda, p, X, Y, offset = rep(0, NROW(X)),
ord = 2, pfixed = FALSE, trace = TRUE)
}

\arguments{
  \item{beta}{Regression parameters}
  \item{lambda}{The scale paramater}
  \item{p}{The shape parameter}
  \item{X}{The design (covariate) matrix.}
  \item{Y}{The response, a survival object.}
  \item{offset}{Offset.}
  \item{ord}{ord = 0 means only lolihood, 1 means score vector as well,
  2 loglihood, score and hessian.}
  \item{pfixed}{Logical, if TRUE the shape parameter is regarded as a known
  constant in the calculations.}
  \item{trace}{If TRUE, more verbose output.}
}
\details{
Note that the function returns log likelihood, score vector
and minus hessian, i.e. the observed information.
The model is
\deqn{h(t; p, \lambda, \beta, z) = \lambda (t / \lambda)^{(p-1)}\exp{-(
  t / \lambda)^p} \exp({z\beta})}
This is in correspondence with \code{\link{dgamma}}.
}
\value{
A list with components
  \item{f}{The log likelihood. Present if \code{ord >= 0}}
  \item{fp}{The score vector. Present if \code{ord >= 1}}
  \item{fpp}{The negative of the hessian. Present if \code{ord >= 2}}
}
%\references{ }
\author{Göran Broström}
%\note{}

\seealso{\code{\link{weibreg}}}
%\examples{
%}
\keyword{survival}
\keyword{distribution}
