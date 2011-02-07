\name{toBinary}
\alias{toBinary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transforms a "survival" data frame into a data frame suitable for
binary (logistic) regression}
\description{
The result of the transformation can be used to do survival analysis via
logistic regression. If the \code{cloglog} link is used, this
corresponds to a discrete time analogue to Cox's proportional hazards model.
}
\usage{
toBinary(dat, surv = c("enter", "exit", "event"),
strats, max.survs = NROW(dat))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame with three variables representing the survival
  response. The default is that they are named \code{enter},
  \code{exit}, and \code{event}} 
  \item{surv}{A character string with the names of the three variables
    representing survival.} 
  \item{strats}{An eventual stratification variable.}
  \item{max.survs}{Maximal numger of survivors per risk set. If set to a
  (small) number, survivors are sampled from the risk sets.}
}
\details{
  toBinary calls \code{risksets} in the \code{eha} package. 
}
\value{
 Returns a data frame expanded risk set by risk set. The three "survival
 variables" are replaced by a variable named \code{event} (which
 overwrites an eventual variable by that name in the input). Two more
 variables are created, \code{riskset} and \code{orig.row}.
  \item{event}{Indicates an event in the corresponding risk set.}
  \item{riskset}{Factor (with levels 1, 2, \ldots) indicating risk
  set.}
  \item{risktime}{The 'risktime' (age) in the corresponding riskset.}
  \item{orig.row}{The row number for this item in the original data frame.}
}
\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}
\note{ 
The survival variables must be three. If you only have \emph{exit} and
\emph{event}, create a third containing all zeros. 
}
\seealso{\code{\link[eha]{coxreg}}, \code{\link[stats]{glm}}.}
\examples{
enter <- rep(0, 4)
exit <- 1:4
event <- rep(1, 4)
z <- rep(c(-1, 1), 2)
dat <- data.frame(enter, exit, event, z)
binDat <- toBinary(dat)
dat
binDat
coxreg(Surv(enter, exit, event) ~ z, method = "ml", data = dat)
## Same as:
summary(glm(event ~ z + riskset, data = binDat, family = binomial(link = cloglog)))
}
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
