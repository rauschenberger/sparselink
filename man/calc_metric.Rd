% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{calc_metric}
\alias{calc_metric}
\title{Calculate deviance}
\usage{
calc_metric(y, y_hat, family)
}
\arguments{
\item{y}{response}

\item{y_hat}{predictor}

\item{family}{character \code{"gaussian"} or \code{"binomial"}}
}
\value{
Returns the deviance (scalar).
}
\description{
Calculates Gaussian deviance (mean-squared error) and binomial deviance.
}
\examples{
n <- 100
family <- "gaussian"
y <- stats::rnorm(n=n)
y_hat <- stats::rnorm(n=n)
calc_metric(y=y,y_hat=y_hat,family=family)

family <- "binomial"
y <- stats::rbinom(n=n,size=1,prob=0.5)
y_hat <- stats::runif(n=n)
calc_metric(y=y,y_hat=y_hat,family=family)

}
\keyword{internal}
