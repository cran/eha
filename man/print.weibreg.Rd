\name{print.weibreg}
\alias{print.weibreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prints weibreg objects}
\description{
The hazard, the cumulative hazard, the density, and the survivor
baseline functions are plotted. 
}
\usage{
print.weibreg(x, digits = max(options()$digits - 4, 3), ...)
}

\arguments{
  \item{x}{A \code{weibreg} object}
  \item{digits}{Precision in printing}
  \item{\dots}{Not used.}
}

\value{
No value is returned.
}
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}
\note{Doesn't work for threeway or higher order interactions. Use
  \code{\link[survival]{print.coxph}} in that case.} 

% ~Make other sections like Warning with \section{Warning }{....} ~ 
 \seealso{\code{\link{weibreg}}, \code{\link[survival]{print.coxph}}}
%\examples{

%}
\keyword{survival}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
