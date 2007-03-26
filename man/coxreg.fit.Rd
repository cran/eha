%\encoding{UTF-8}
\name{coxreg.fit}
\alias{coxreg.fit}

\title{Cox regression}
\description{
Called by \code{\link{coxreg}}, but a user can call it directly.
}
\usage{
coxreg.fit(X, Y, rs, strats, offset, init, max.survs,
method = "breslow", boot = FALSE, calc.hazards = TRUE,
calc.martres = TRUE, control, verbose = TRUE)
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
  \item{method}{Either "efron" (default) or "breslow".}
  \item{boot}{Number of bootstrap replicates. Defaults to FALSE, no
    bootstrapping.}
  \item{calc.hazards}{Should estimates of baseline hazards be
    calculated?}
  \item{calc.martres}{Should martingale residuals be calculated?}
  \item{control}{See \code{\link{coxreg}}}
  \item{verbose}{Should Warnings about convergence be printed?}
}
\details{
\code{rs} is dangerous to use when NA's are present. 
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
  \item{bootstrap}{The bootstrap replicates, if requested on input.}
  \item{conver}{\code{TRUE} if convergence.}
  \item{f.conver}{TRUE if variables converged.}
  \item{fail}{\code{TRUE} if failure.}
  \item{iter}{Number of performed iterations.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
\note{It is the user's responsibility to check that indata is sane.}

\seealso{\code{\link{coxreg}}, \code{\link{risksets}}}
\examples{
 X <- as.matrix(data.frame(
                x=     c(0, 2,1,4,1,0,3),
                sex=   c(1, 0,0,0,1,1,1)))
 time <- c(1,2,3,4,5,6,7)
 status <- c(1,1,1,0,1,1,0)
 stratum <- rep(1, length(time))

 coxreg.fit(X, Surv(time, status), strats = stratum, max.survs = 6,
     control = list(eps=1.e-4, maxiter = 10, trace = FALSE))
}

\keyword{survival}
\keyword{regression}
