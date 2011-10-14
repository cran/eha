\name{ltx}
\alias{ltx}
\title{
LaTeX printing of coxreg results.
}
\description{
This function prints the LaTeX code of the results of a fit from
\code{\link{coxreg}}, similar to what \code{xtable} does for fits from
other functions.
}
\usage{
ltx(x, caption = NULL, label = NULL, dr = NULL, digits = max(options()$digits - 4, 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The output from a call to \code{coxreg}
}
  \item{caption}{
A suitable caption for tha table.
}
  \item{label}{
A label used in the LaTeX code.
}
  \item{dr}{
Output from a \code{drop1} call with input the fit from a call to \code{coxreg}.
}
  \item{digits}{
Number of digits to be printed.
}
  \item{\dots}{
Not used.
}
}
\details{
The result is a printout which is (much) nicer than the standard printed
output from \code{glm} and friends,
}
\value{
LaTeX code version of the results from a run with \code{\link{coxreg}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Göran Broström.
}
\note{
There is no method in \code{xtable} for \code{coxreg}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{xtable}, \code{\link{coxreg}}
}
\examples{
data(oldmort)
fit <- coxreg(Surv(enter, exit, event) ~ civ + sex, data = oldmort)
dr <- drop1(fit, test = "Chisq")
ltx(fit, dr = dr, caption = "A test example.", label = "tab:test1") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{printing}

