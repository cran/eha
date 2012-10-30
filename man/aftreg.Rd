\name{aftreg}
\alias{aftreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accelerated Failure Time Regression}
\description{The accelerated failure time model with parametric baseline
  hazard(s).  Allows for stratification with different
  scale and shape in each stratum, and left truncated and right censored data.
}
\usage{
aftreg(formula = formula(data), data = parent.frame(),
na.action = getOption("na.action"), dist = "weibull", init, shape = 0,
id, param = c("default", "survreg", "canonical"),
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
    alternatives "gompertz", "ev", "loglogistic" and
    "lognormal". A special 
    case like the \code{exponential} can be obtained by choosing "weibull"
    in combination with \code{shape = 1}.}
  \item{init}{vector of initial values of the iteration.
    Default initial value is zero for all variables.}
  \item{shape}{If positive, the shape parameter is fixed at that value.
    If zero or negative, the shape parameter is estimated. Stratitfication
    is not meaningful if shape is fixed.}
  \item{id}{If there are more than one spell per individual, it is
    essential to keep spells together by the id argument. This allows for
    time-varying covariates.}
  \item{param}{Which parametrization should be used? The \code{default}
    uses the parametrization given in the vignette, while the
    \code{survreg} uses the same as in the \code{\link[survival]{survreg}}
      function. The \code{canonical} separates the
      parameters according to AFT or PH. Only(?) useful for the
      \code{Gompertz} distribution.}
    \item{control}{a list with components \code{eps} (convergence
      criterion), \code{maxiter} (maximum number of iterations), and
      \code{trace} (logical, debug output if \code{TRUE}).  You can
      change any component without mention the other(s).}
    \item{singular.ok}{Not used.}
    \item{model}{Not used.}
    \item{x}{Return the design matrix in the model object?}
    \item{y}{Return the response in the model object?}
    \item{center}{Deprecated. No centering results are reported.}
  }

\details{
The parameterization is different
from the one used by \code{\link[survival]{survreg}}, when \code{param =
  "default"}.
Then the model is
\deqn{S(t; a, b, \beta, z) = S_0((t / \exp(b -
  z\beta))^{\exp(a)})}{%
  S(t; a, b, beta, z) =  S_0((t/exp(b - z beta))^exp(a))}
where \eqn{S_0} is some standardized survivor function.
The base-line parameters \eqn{a} and \eqn{b} are log shape and log
scale, respectively. This is for the \code{default}
parametrization. With the \code{survreg} parametrization, some signs are
changed: \deqn{b - z beta} is changed to \deqn{b + z beta}, and
\deqn{exp(a)} to \deqn{exp(-a)}. For \code{canonical}, in addition,
shape is divided by scale, which implies that the new scale indexes an
AFT family for each fixed shape, while the new shape (like the old)
indexes a PH family for each fixed scale, but only for the Gompertz
distribution.    
}

\value{
A list of class \code{c("aftreg", "coxreg")} with components
  \item{coefficients}{Fitted parameter estimates.}
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
  \item{param}{The parametrization.
  }
}

\author{\enc{Göran Broström}{George Bridgewater}}


%\section{Warning}{The print method \code{\link{print.aftreg}} doesn't work if
%threeway or higher order interactions are present. Use
%\code{\link[survival]{print.coxph}} in that case.
%}

\seealso{\code{\link{coxreg}}, \code{\link{phreg}}, \code{link[survival]{survreg}}
  }

\examples{
data(mort)
aftreg(Surv(enter, exit, event) ~ ses, data = mort)
}
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
