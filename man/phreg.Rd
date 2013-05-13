\name{phreg}
\alias{phreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric Proportional Hazards Regression}
\description{Proportional hazards model with parametric baseline
  hazard(s).  Allows for stratification with different 
  scale and shape in each stratum, and left truncated and right censored data.
}
\usage{
phreg(formula = formula(data), data = parent.frame(),
na.action = getOption("na.action"), dist = "weibull", cuts = NULL,
init, shape = 0, param = c("canonical", "rate"),
control = list(eps = 1e-08, maxiter = 20, trace = FALSE),
singular.ok = TRUE, model = FALSE, x = FALSE, y = TRUE, center = NULL)
}

\arguments{
  \item{formula}{a formula object, with the response on the left of
    a ~ operator, and the terms on the right.  The response must be
    a survival object as returned by the Surv function.}
  \item{data}{a data.frame in which to interpret the variables named
    in the formula.}
  \item{na.action}{ a missing-data filter function, applied to the
    model.frame, after any subset argument has been used.
    Default is \code{options()$na.action}.}
\item{dist}{Which distribution? Default is "weibull", with the
  alternatives "ev" (Extreme value), "gompertz", "pch" (piecewise
  constant hazards function), "loglogistic" and "lognormal". A special 
  case like the \code{exponential} can be obtained by choosing "weibull"
  in combination with \code{shape = 1}.}
\item{cuts}{Only used with \code{dist = "pch"}. Specifies the points in
  time where the hazard function jumps.}
  \item{init}{vector of initial values of the iteration.
    Default initial value is zero for all variables.}
  \item{shape}{If positive, the shape parameter is fixed at that value
    (in each stratum).  If zero or negative, the shape parameter is
    estimated.  If more than one stratum is present in data, each stratum
    gets its own estimate.}
  \item{param}{Applies only to the Gompertz distribution: "canonical" is
  defined in the description of the \code{\link{Gompertz}} distribution;
  "rate" transforms \code{scale} to 1/log(scale), giving the same
  parametrization as in Stata and SAS. The latter thus allows for a
  negative rate, or a "cure" (Gompertz) model. The default is
  "canonical"; if this results in extremely large scale and/or shape
  estimates, consider trying "rate".}
    \item{control}{a list with components \code{eps} (convergence
    criterion), \code{maxiter} (maximum number of iterations), and
    \code{silent} (logical, controlling amount of output).  You can
    change any component without mention the other(s).}
  \item{singular.ok}{Not used.}
  \item{model}{Not used.}
  \item{x}{Return the design matrix in the model object?}
  \item{y}{Return the response in the model object?}
  \item{center}{Deprecated (not used); results are now reported as is, without
  centering. You may center your covariates manually, if you want centered
  results. Centering affects only 'scale' and interaction terms.} 
}

\details{
The parameterization is the same as in \code{\link{coxreg}}
and \code{\link[survival]{coxph}}, but different
from the one used by \code{\link[survival]{survreg}} (which is not a
proportional hazards modelling function).
The model is
\deqn{S(t; a, b, \beta, z) = S_0((t/b)^a)^{\exp((z-mean(z))\beta)}}{%
  S(t; a, b, beta, z) =  S0((t/b)^a)^exp((z - mean(z)) beta)}
where S0 is some standardized survivor function.
}

\value{
A list of class \code{c("phreg", "coxreg")} with components
\item{coefficients}{Fitted parameter estimates.}
\item{cuts}{Cut points for the "pch" distribution. \code{NULL}
  otherwise.}
\item{hazards}{The estimated constant levels in the case of the "pch"
  distribution. \code{NULL} otherwise.}
  \item{var}{Covariance matrix of the estimates.}
  \item{loglik}{Vector of length two; first component is the value at
    the initial parameter values, the second componet is the maximized value.}
  \item{score}{The score test statistic (at the initial value).}
  \item{linear.predictors}{The estimated linear predictors.}
  \item{means}{Means of the columns of the design matrix.}
  \item{w.means}{Weighted (against exposure time) means of covariates;
    weighted relative frequencies of levels of factors.}
  \item{n}{Number of spells in indata (possibly after removal of cases
    with NA's).}
  \item{events}{Number of events in data.}
  \item{terms}{Used by extractor functions.}
  \item{assign}{Used by extractor functions.}
%  \item{wald.test}{The Wald test statistic (at the initial value).}
  \item{y}{The Surv vector.}
  \item{isF}{Logical vector indicating the covariates that are factors.}
  \item{covars}{The covariates.}
  \item{ttr}{Total Time at Risk.}
  \item{levels}{List of levels of factors.}
  \item{formula}{The calling formula.}
  \item{call}{The call.}
  \item{method}{The method.}
  \item{convergence}{Did the optimization converge?}
  \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
  \item{pfixed}{TRUE if shape was fixed in the estimation.}
}

\author{\enc{Göran Broström}{George Bridgewater}}
\note{For some of the baseline distributions, e.g., the lognormal
  distribution, the family of distributions is not closed under
  proportional hazards.}


\section{Warning}{Note that covariates are internally centered by this
  function,  
and this is not corrected for in the output. This affects the estimate
of \eqn{\log(scale)}{log(scale)}, but nothing else. It is possible to
  choose not to center, though. }

\seealso{\code{\link{coxreg}}, \code{\link{check.dist}}, \code{link{aftreg}}.}

\examples{
data(mort)
fit <- phreg(Surv(enter, exit, event) ~ ses, data = mort)
fit
plot(fit)
fit.cr <- coxreg(Surv(enter, exit, event) ~ ses, data = mort)
check.dist(fit.cr, fit)
}
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
