\name{Coxreg}
\alias{Coxreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cox regression}
\description{
  Performs Cox regression with some special attractions, especially
  \emph{sampling of risksets} and \emph{the weird bootstrap}.
}
\usage{
Coxreg(formula = formula(data), data = parent.frame(), weights, subset,
na.action = getOption("na.action"), init = NULL,
method = c("efron", "breslow"),
baseline = FALSE,
control = list(eps = 1e-08, maxiter = 25, trace = FALSE),
singular.ok = TRUE, model = FALSE,
x = FALSE, y = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a formula object, with the response on the left of
    a ~ operator, and the terms on the right. The response must be
    a survival object as returned by the Surv function.}
  \item{data}{a data.frame in which to interpret the variables named
    in the formula.}
  \item{weights}{Case weights; time-fixed or time-varying.}
  \item{subset}{An optional vector specifying a subset of observations
    to be used in the fitting process.} 
  \item{na.action}{ a missing-data filter function, applied to the
    model.frame, after any subset argument has been used.
    Default is \code{options()$na.action}.}
  \item{init}{vector of initial values of the iteration.
    Default initial value is zero for all variables.}
  \item{method}{Method of treating ties, "efron" (default), "breslow".}
  \item{baseline}{Should the baseline hazards be estimated? Default is No.}
  \item{control}{a list with components \code{eps} (convergence
    criterion), \code{maxiter} (maximum number of iterations), and
    \code{silent} (logical, controlling amount of output). You can
    change any component without mention the other(s).}
  \item{singular.ok}{Not used}
  \item{model}{Not used}
  \item{x}{Return the design matrix in the model object?}
  \item{y}{return the response in the model object?}
}
\details{
The default method, \code{efron}, and the alternative, \code{breslow},
  are both the same as in \code{\link[survival]{coxph}} in package
  \code{survival}. 
}
\value{
A list of class \code{c("coxreg", "coxph")} with components
  \item{coefficients}{Fitted parameter estimates.}
  \item{var}{Covariance matrix of the estimates.}
  \item{loglik}{Vector of length two; first component is the value at
    the initial parameter values, the second componet is the maximized value.}
  \item{score}{The score test statistic (at the initial value).}
  \item{linear.predictors}{The estimated linear predictors.}
  \item{residuals}{The martingale residuals.}
  \item{hazard}{The estimated baseline hazard, calculated at the means
    of the covariates (rather, columns of the design matrix). Is a list,
    with one component per stratum. Each 
    component is a matrix with two columns, the first contains
    risktimes, the second the corresponding hazard atom.}
  \item{means}{Means of the columns of the design matrix.}
  \item{w.means}{Weighted (against exposure time) means of covariates;
    weighted relative frequencies of levels of factors.}
  \item{n}{Number of spells in indata (possibly after removal of cases
    with NA's).}
  \item{events}{Number of events in data.}
  \item{terms}{Used by extractor functions.}
  \item{assign}{Used by extractor functions.}
  \item{wald.test}{The Wald test statistic (at the initial value).}
  \item{y}{The Surv vector.}
  \item{isF}{Logical vector indicating the covariates that are factors.}
  \item{covars}{The covariates.}
  \item{ttr}{Total Time at Risk.}
  \item{levels}{List of levels of factors.}
  \item{formula}{The calling formula.}
  \item{bootstrap}{The (matrix of) bootstrap replicates, if requested on
    input. It is up to the user to do whatever desirable with this sample.}
  \item{boot.sd}{The estimated standard errors of the bootstrap replicates.}
  \item{call}{The call.}
  \item{method}{The method.}
  \item{convergence}{Did the optimization converge?}
  \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
}
\references{Broström, G. and Lindkvist, M. (2008). Partial partial
  likelihood. Communications in Statistics: Simulation and Computation
  37:4, 679-686.}

\author{\enc{Göran Broström}{George Bridgewater}}

\note{This function starts by creating risksets, if no riskset is
  supplied via \code{rs}, with the aid of
  \code{\link{risksets}}. Supplying output from \code{risksets} via \code{rs}
  fails if there are any NA's in the data! Note also
  that it depends on stratification, so \code{rs} contains information
  about stratification. Giving another strata variable in the formula is
  an error. The same is ok, for instance to supply stratum
  interactions.}
\section{Warning}{The use of \code{rs} is dangerous, see note. It
  can however speed up computing time considerably for huge data sets.} 

\seealso{\code{\link[survival]{coxph}}, \code{\link{risksets}}}
\examples{

 dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
                status=c(1,1,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
 Coxreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
 }
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
