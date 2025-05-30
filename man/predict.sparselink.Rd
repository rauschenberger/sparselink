% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{predict.sparselink}
\alias{predict.sparselink}
\title{Out-of-sample Predictions}
\usage{
\method{predict}{sparselink}(object, newx, weight = NULL, ...)
}
\arguments{
\item{object}{object of class \code{"sparselink"}
(generated by function \link{sparselink})}

\item{newx}{features:
matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
for multi-task learning;
list of \eqn{q} matrices
with \eqn{n_k} rows (samples) and \eqn{p} columns (variables)
for transfer learning, for each \eqn{k} in \eqn{1,\ldots,q}}

\item{weight}{hyperparameters for scaling external and internal weights:
numeric vector of length 2,
with the first entry for the external weights
(prior coefficients from source data),
and the second entry for the internal weights
(prior coefficients from target data),
selected values must be among the candidate values,
default: \code{NULL} (using cross-validated weights)}

\item{...}{(not applicable)}
}
\value{
Returns predicted values or predicted probabilities.
The output is a list of \eqn{q} column vectors of length \eqn{n_k}
for \eqn{k} in \eqn{1,\ldots,q}.
Each vector corresponds to one target (multi-task learning)
or one dataset (transfer learning).
}
\description{
Predicts outcomes with a multi-task or transfer learning regression model.
}
\examples{
family <- "gaussian"
type <- "multiple" # try "multiple" or "transfer"
if(type=="multiple"){
 data <- sim_data_multi(family=family)
} else if(type=="transfer"){
 data <- sim_data_trans(family=family)
}
\dontshow{data <- sim_data_multi(family=family,n0=10,p=3)}
object <- sparselink(x=data$X_train,y=data$y_train,family=family)
y_hat <- predict(object=object,newx=data$X_test)

}
\references{
\href{https://orcid.org/0000-0001-6498-4801}{Armin Rauschenberger},
\href{https://orcid.org/0000-0003-3443-0298}{Petr N. Nazarov}, and
\href{https://orcid.org/0000-0003-3977-7469}{Enrico Glaab}
(2025).
"Estimating sparse regression models in multi-task learning and transfer learning through adaptive penalisation".
\emph{Under revision}. \url{https://hdl.handle.net/10993/63425}
}
\seealso{
Use \code{\link{sparselink}} to fit the model
and \code{\link[=coef.sparselink]{coef}} to extract coefficients.
}
