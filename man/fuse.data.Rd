% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{fuse.data}
\alias{fuse.data}
\title{Data fusion}
\usage{
fuse.data(x, y = NULL, foldid = NULL)
}
\arguments{
\item{x}{list of q matrices, with n_1,...,n_q rows and with p columns}

\item{y}{list of q vectors, of length n_1,...,n_q, or NULL (default)}

\item{foldid}{list of q vectors, of length n_1,...n_q, or NULL (default)}
}
\description{
Data fusion
}
\examples{
data <- sim.data.transfer()
sapply(X=data$y_train,FUN=length)
sapply(X=data$X_train,FUN=dim)
fuse <- fuse.data(x=data$X_train,y=data$y_train)
length(fuse$y)
dim(fuse$x)
table(fuse$index)

}
