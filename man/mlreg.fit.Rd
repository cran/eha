\name{mlreg.fit}
\alias{mlreg.fit}

\title{ML proportional hazards regression}
\description{
Called by \code{\link{mlreg}}, but a user can call it directly.
}
\usage{
mlreg.fit(X, Y, rs, strats, offset, init, max.survs,
method = "ML", boot = FALSE, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{The design matrix.}
  \item{Y}{The survival object.}
  \item{rs}{The risk set composition. If absent, calculated.}
  \item{strats}{The stratum variable. Can be absent.}
  \item{offset}{Offset. Can be absent.}
  \item{init}{Start values. If absent, equal to zero.}
  \item{max.survs}{Sampling of risk sets? If so, gives the maximum
    number of survivors in each risk set.}
  \item{method}{Either "ML" (default) or "MPPL".}
    \item{boot}{No. of bootstrap replicates. Defaults to FALSE, i.e., no
    bootstrapping.} 
  \item{control}{See \code{\link{coxreg}}}
}
\details{
  See \code{\link{mlreg}} for details.
}
\value{
  A list with components
  \item{coefficients}{Estimated regression parameters.}
  \item{var}{Covariance matrix of estimated coefficients.}
  \item{loglik}{First component is value at \code{init}, second at maximum.}
  \item{score}{Score test statistic, at initial value.}
  \item{linear.predictors}{Linear predictors.}
  \item{residuals}{Martingale residuals.}
  \item{hazard}{Estimated baseline hazard. At value zero of
    'design' variables.}
  \item{means}{Means of the columns of the design matrix.}
  \item{bootstrap}{The bootstrap sample, if requested on input.}
  \item{conver}{\code{TRUE} if convergence.}
  \item{fail}{\code{TRUE} if failure.}
  \item{iter}{Number of performed iterations.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
\note{\code{rs} is dangerous to use when NA's are present. 
It is the user's responsibility to ensure that indata is sane.}

\seealso{\code{\link{coxreg}}, \code{\link{risksets}}}
\examples{
 X <- as.matrix(data.frame(
                x=     c(0, 2,1,4,1,0,3),
                sex=   c(1, 0,0,0,1,1,1)))
 time <- c(1,2,3,4,5,6,7)
 status <- c(1,1,1,0,1,1,0)
 stratum <- rep(1, length(time))

 mlreg.fit(X, Surv(time, status), strats = stratum, max.survs = 6,
     control = list(eps=1.e-4, maxiter = 10, trace = TRUE))
}

\keyword{survival}
\keyword{regression}
