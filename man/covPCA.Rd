\name{covPCA}
\alias{covPCAproj}
\alias{covPCAgrid}
\title{ Robust Covariance Matrix Estimation }
\description{
  computes the robust covariance matrix using the \code{PCAgrid} and 
  \code{PCAproj} functions.
}
\usage{
covPCAproj(x, control, ...)
covPCAgrid(x, control, ...)
}
\arguments{
  \item{x}{ a numeric matrix or data frame which provides the data. }
  \item{control}{ a list whose elements must be the same as (or a subset of) 
  the parameters of the appropriate PCA function (\code{\link{PCAgrid}} or 
  \code{\link{PCAproj}}). If the control object is supplied, other parameters 
  supposed to be passed to this function as part of the \code{...} parameter 
  are overridden. }
  \item{\dots}{ additional arguments passed to the function 
  \code{\link{PCAproj}} or \code{\link{PCAgrid}} respectively.}
}
\details{
  The functions \code{covPCAproj} and \code{covPCAgrid} use the functions
  \code{\link{PCAproj}} and \code{\link{PCAgrid}} respectively to estimate 
  the covariance matrix of the data matrix \code{x}.
}
\value{
  \item{cov}{ the actual covariance matrix estimated from \code{x}}
  \item{center}{ the center of the data \code{x} that was substracted from them
  before the PCA algorithms were run.}
  \item{method}{ a string describing the method that was used to calculate the
  covariance matrix estimation}
}
\references{ 
  C. Croux, P. Filzmoser, M. Oliveira (2004)
  \emph{Projection-pursuit Estimators for Robust Principal Component Analysis, 
        Technical Report TS-04-4, Vienna University of Technology, Austria}
}
\author{ Heinrich Fritz, Peter Filzmoser <\email{P.Filzmoser@tuwien.ac.at}> }
\seealso{ \code{\link{PCAgrid}}, \code{\link{ScaleAdvR}}, \code{\link{princomp}} }
\examples{
  # multivariate data with outliers
  x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
             rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
  covPCAproj(x)
  # compare with classical covariance matrix:
  cov(x)
}
\keyword{ robust }
\keyword{ multivariate }
