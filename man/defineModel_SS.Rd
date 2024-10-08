% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defineModel_SS.R
\name{defineModel_SS}
\alias{defineModel_SS}
\title{defineModel_SS}
\usage{
defineModel_SS(inclPhenology = TRUE, ListLen = NULL, inclStateRE = FALSE)
}
\arguments{
\item{inclPhenology}{should the model account for seasonal variation?}

\item{ListLen}{should list length be modelled as continuous (\code{cont}), categorical (\code{cat}) or ignored (\code{NULL})}

\item{inclStateRE}{should there be a site-level random effect in the state model?}
}
\value{
a set of code
}
\description{
defineModel_SS
}
\details{
Defines the occupancy model in Nimble for one species. Currently this implements a simple occupancy model.
This is very close to the random walk formulation used in sparta
There are some minor differences: priors are generally more constrained and \code{alpha.0} is simplified (\code{alpha.p} in sparta)
}
