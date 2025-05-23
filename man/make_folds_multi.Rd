% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{make_folds_multi}
\alias{make_folds_multi}
\alias{make_folds_trans}
\title{Create folds for multi-task and transfer learning}
\usage{
make_folds_multi(y, family, nfolds = 10)

make_folds_trans(y, family, nfolds = 10)
}
\arguments{
\item{y}{multi-task learning:
y matrix with \eqn{n} rows (samples) and \eqn{q} columns (outcomes)
transfer learning:
list of \eqn{q} numeric vectors}

\item{family}{character \code{"gaussian"} or \code{"binomial"}}

\item{nfolds}{integer between 2 and \eqn{n}}
}
\description{
Create folds for multi-task and transfer learning
}
\examples{
#--- multi-task learning ---
family <- "binomial"
y <- sim_data_multi(family=family)$y_train
fold <- make_folds_multi(y=y,family=family)

#--- transfer learning ---
family <- "binomial"
y <- sim_data_trans(family=family)$y_train
fold <- make_folds_trans(y,family=family)

}
\keyword{internal}
