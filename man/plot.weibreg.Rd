\name{plot.weibreg}
\alias{plot.weibreg}

\title{Plots output from a Weibull regression}
\description{
Just a simple plot of the hazard functions for each stratum.
}
\usage{
plot.weibreg(x, new.data = rep(0, length(x$means)), ...)
}

\arguments{
  \item{x}{A \code{weibreg} object}
  \item{new.data}{At which covariate values? This doesn't work for the moment.}
  \item{...}{Who knows?}
}
\warning{\code{new.data} does not work as intended.}
\details{
The plot is drawn at the mean values of the covariates.
}
\value{
No return
}

\author{Göran Broström}

\seealso{\code{\link{weibreg}}, ~~~ }
\examples{
## The function is currently defined as
function (x, new.data = rep(0, length(x$means))) 
{
    if (!inherits(x, "weibreg")) 
        stop("Works only with 'weibreg' objects.")
    ncov <- length(x$means)
    ns <- x$n.strata
    lambda <- exp(x$coefficients[ncov + (1:ns) * 2 - 1])
    p <- exp(x$coefficients[ncov + (1:ns) * 2])
    xlim <- c(min(x$y[, 1]), max(x$y[, 2]))
    npts <- 199
    xx <- seq(xlim[1], xlim[2], length = npts)
    if (xx[1] <= 0) 
        xx[1] <- 0.001
    haz <- matrix(ncol = npts, nrow = ns)
    for (i in 1:ns) {
        tl <- xx/lambda[i]
        haz[i, ] <- (p[i]/lambda[i]) * tl^(p[i] - 1)
    }
    ylim <- c(0, max(haz))
    cat("length(xx) = ", length(xx), "\n")
    cat("dim(haz) = ", dim(haz), "\n")
    plot(xx, haz[1, ], type = "l", xlim = xlim, ylim = ylim)
    if (ns > 1) {
        for (i in 2:ns) {
            lines(xx, haz[i, ], type = "l", lty = i)
        }
    }
  }
}
\keyword{dplot}% at least one, from doc/KEYWORDS
\keyword{survival}% __ONLY ONE__ keyword per line
