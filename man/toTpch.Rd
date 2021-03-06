% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toTpch.R
\name{toTpch}
\alias{toTpch}
\title{Transform survival data to tabular form}
\usage{
toTpch(formula, data, cuts, enter = "enter", exit = "exit",
event = "event", episode = "age")
}
\arguments{
\item{formula}{A model formula.}

\item{data}{A data frame with survival data.}

\item{cuts}{An ordered, non-negative vector of time points at which a hazard function changes value.
Note that data are left truncated at cuts[1] (the smallest value) and right censored at c[n], where 
n  is the length of cuts and cuts[n] == max(cuts).}

\item{enter}{Character string with the name of the variable representing left truncation values.}

\item{exit}{Character string with the name of the event/censoring time variable.}

\item{event}{Character string with the name of the event indicator variable.}

\item{episode}{Character string with the name of the output variable of the grouped time (a factor variable).}
}
\value{
A data frame with exposure time and number of events aggregated by time intervals and covariates.
 If all covariates are factors,this usually results in a huge reduction of the size of thedata frame, 
 but otherwise the size of the output may be larger than the size of the input data frame
}
\description{
Transform a "survival data frame" to tabular form aggregating number of events 
and exposure time by time intervals and covariates.
}
\details{
If \code{cuts} is missing, nothing is done. Internally, this function first calls 
 \code{survival::survSplit} and then \code{stats::aggregate}.
}
\note{
Episodes, or parts of episodes, outside \code{min(cuts), max(cuts)} are cut off. 
 With continuous covariates, consider rounding them so that the number of distinct oberved values is not too large.
}
\author{
Göran Broström
}
\keyword{manip}
