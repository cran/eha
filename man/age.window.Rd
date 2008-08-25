\name{age.window}
\alias{age.window}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Age cut of survival data}
\description{For a given age interval, each spell is cut to fit into the
  given age interval.
}
\usage{
age.window(dat, window, surv=c("enter", "exit", "event"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{Input data frame. Must contain survival data.}
  \item{window}{Vector of length two; the age interval.}
  \item{surv}{Vector of length three giving the names of the central
    variables in 'dat'.}
}
\details{
The \code{window} must be in the order \code{(begin, end)}
}
\value{
  A data frame of the same form as the input data frame, but 'cut' as
  desired. Intervals exceeding \code{window[2]} will be given
  \code{event = 0}

}
%\references{ ~put references to the literature/web site here ~ }
\author{Göran Broström}
%\note{ ~~further notes~~ }

\seealso{\code{\link{cal.window}}, \code{\link{coxreg}},
\code{\link{aftreg}}}

\examples{
dat <- data.frame(enter = 0, exit = 5.731, event = 1, x = 2)
window <- c(2, 5.3)
dat.trim <- age.window(dat, window)  
}
\keyword{survival}% at least one, from doc/KEYWORDS
