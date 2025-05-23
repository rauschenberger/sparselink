% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{count_vector}
\alias{count_vector}
\alias{count_matrix}
\title{Metrics for sign detection}
\usage{
count_vector(truth, estim)

count_matrix(truth, estim)
}
\arguments{
\item{truth}{(i) vector of length \eqn{p} or
(ii) \eqn{n \times p} matrix with entries in -1, 0, 1}

\item{estim}{(i) vector of length \eqn{p} or
(ii) \eqn{n \times p} matrix with entries -1, 0, 1}
}
\description{
Calculates sensitivity, specificity and precision for ternary data
(with -1 for negative effect, 0 for no effect, 1 for positive effect).
}
\examples{
truth <- sample(x=c(-1,0,1),size=20,replace=TRUE)
estim <- sample(x=c(-1,0,1),size=20,replace=TRUE)
table(truth,estim)
count_vector(truth,estim)

}
\keyword{internal}
