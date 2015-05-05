\name{print.risksets}
\alias{print.risksets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints a summary of the content of a set of risk sets.
}
\description{
Given the output from \code{risksets}, summary statistics are given for it.
}
\usage{
\method{print}{risksets}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class 'risksets'.
}
  \item{...}{Not used for the moment.}
}

\value{
No value is returned; the function prints summary statistics of risk sets.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Göran Broström
}
\note{
There is no \code{summary.risksets} yet. On the TODO list.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{risksets}}
%% ~~objects to See Also as \code{\link{help}}, ~~~

\examples{
rs <- with(mort, risksets(Surv(enter, exit, event)))
print(rs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary}

