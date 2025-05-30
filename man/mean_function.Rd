% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mean_function}
\alias{mean_function}
\title{Mean function}
\usage{
mean_function(eta, family)
}
\arguments{
\item{eta}{numeric vector}

\item{family}{character \code{"gaussian"} or \code{"binomial"}}
}
\value{
Returns a numeric vector of the transformed values.
}
\description{
Applies the mean function (inverse link function).
}
\examples{
family <- "binomial"
eta <- seq(from=-3,to=3,length.out=100)
mu <- mean_function(eta=eta,family=family)
graphics::plot(x=eta,y=mu,type="l",main=family)
graphics::abline(v=0,lty=2)
h <- ifelse(family=="binomial",0.5,0)
graphics::abline(h=h,lty=2)

}
\keyword{internal}
