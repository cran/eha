\encoding{UTF-8}
\name{print.coxreg}
\alias{print.coxreg}

\title{Prints coxreg objects}
\description{
More "pretty-printing" than \code{print.coxph}, which is a fall-back for
'difficult' objects.
}
\usage{
print.coxreg(x, digits = max(options()$digits - 4, 3), ...)
}
\arguments{
  \item{x}{A \code{coxreg} object, typically the result of running
    \code{coxreg}} 
  \item{digits}{Output format.}
  \item{\dots}{Other arguments.}
}
\details{
Doesn't work with three-way and higher interactions, in which case
\code{print.coxph} is used. Prints also output from \code{\link{mlreg}}.
}
\value{
No value is returned.
 }
%\references{ ~put references to the literature/web site here ~ }
\author{\enc{Göran Broström}{George Bridgewater}}
%\note{}


\seealso{\code{\link{coxreg}}, \code{\link[survival]{print.coxph}}}
%\examples{

}
\keyword{survival}% at least one, from doc/KEYWORDS
